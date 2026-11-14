use crate::cdc_soda::now_unix_secs;
use epi_core::{
    parse_flexible_date_to_mmwr, EpiError, MmwrWeek, SourceKind, TimeSeriesDataset, TimeSeriesPoint,
};
use epi_store::{
    encode_pprl_record, link_and_deduplicate_pprl, EncodedPprlRecord, EpiStore, PprlConfig,
    PprlDemographicInput, PprlLinkageSummary,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

/// Legacy & Hospital Public Health Messaging: Pipe-Delimited HL7 Version 2.5.1
/// Electronic Laboratory Reporting (`ORU^R01`) and Emergency Department Syndromic Surveillance (`ADT^A01/A03/A04/A08`)
/// coupled with inline 512-bit `HMAC-SHA256` Bloom-filter `PPRL` patient deduplication.
///
/// References:
/// - HL7 Version 2.5.1 Implementation Guide: Electronic Laboratory Reporting (ELR) to Public Health: https://www.hl7.org/implement/standards/product_brief.cfm?product_id=98
/// - CDC National Syndromic Surveillance Program (NSSP) HL7 v2.5.1 Implementation Guides: https://www.cdc.gov/nssp/php/onboarding-resources/implementation-guides.html
/// - NSSP Syndromic Surveillance Messaging Architecture: https://pmc.ncbi.nlm.nih.gov/articles/PMC6606111/

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NsspSyndromeCategory {
    IliRespiratory,
    RashExanthem,
    GiEnteric,
    Neurologic,
    UnspecifiedOther,
}

impl NsspSyndromeCategory {
    pub fn label(self) -> &'static str {
        match self {
            Self::IliRespiratory => "ILI / Respiratory",
            Self::RashExanthem => "Rash / Exanthem",
            Self::GiEnteric => "GI / Enteric",
            Self::Neurologic => "Neurologic",
            Self::UnspecifiedOther => "Other / Unclassified",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParsedHl7v2Event {
    pub message_control_id: String,
    pub message_type: String,
    pub sending_facility: String,
    pub event_iso_date: String,
    pub loinc_or_chief_code: String,
    pub result_or_complaint_text: String,
    pub is_positive_or_case: bool,
    pub syndrome: NsspSyndromeCategory,
    pub state_abbr: String,
    /// Zero-PII 16-hex-char prefix of the 512-bit HMAC-SHA256 Bloom filter for UI audit.
    pub bloom_fingerprint_prefix: String,
    pub linked_entity_id: usize,
    pub retained_after_pprl: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Hl7v2IngestOutcome {
    pub dataset: TimeSeriesDataset,
    pub total_messages_parsed: usize,
    pub positive_or_syndromic_hits: usize,
    pub deduplicated_incident_cases: usize,
    pub syndrome_counts: Vec<(String, usize)>,
    pub pprl_summary: PprlLinkageSummary,
    pub events: Vec<ParsedHl7v2Event>,
}

fn hl7_ts_to_iso_and_days(ts: &str) -> Option<(String, f64)> {
    let digits: String = ts.chars().take_while(|c| c.is_ascii_digit()).collect();
    if digits.len() < 8 {
        return None;
    }
    let y: i64 = digits[0..4].parse().ok()?;
    let m: i64 = digits[4..6].parse().ok()?;
    let d: i64 = digits[6..8].parse().ok()?;
    let iso = format!("{y:04}-{m:02}-{d:02}");

    let y_adj = y - if m <= 2 { 1 } else { 0 };
    let era = if y_adj >= 0 { y_adj } else { y_adj - 399 } / 400;
    let yoe = y_adj - era * 400;
    let doy = (153 * (m + if m > 2 { -3 } else { 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    let abs_days = (era * 146097 + doe - 719468) as f64;
    Some((iso, abs_days))
}

fn classify_syndrome(text: &str, code: &str) -> NsspSyndromeCategory {
    let u = format!("{} {}", text.to_ascii_uppercase(), code.to_ascii_uppercase());
    if u.contains("RASH")
        || u.contains("MACULOPAPULAR")
        || u.contains("MEASLES")
        || u.contains("B05")
        || u.contains("14189004")
    {
        NsspSyndromeCategory::RashExanthem
    } else if u.contains("COUGH")
        || u.contains("FEVER")
        || u.contains("DYSPNEA")
        || u.contains("INFLUENZA")
        || u.contains("PERTUSSIS")
        || u.contains("J09")
        || u.contains("J11")
    {
        NsspSyndromeCategory::IliRespiratory
    } else if u.contains("VOMIT") || u.contains("DIARRHEA") || u.contains("GASTRO") || u.contains("A08") {
        NsspSyndromeCategory::GiEnteric
    } else if u.contains("MENING") || u.contains("ENCEPH") || u.contains("SEIZURE") {
        NsspSyndromeCategory::Neurologic
    } else {
        NsspSyndromeCategory::UnspecifiedOther
    }
}

fn is_elr_result_positive(obx5: &str, obx8_flag: &str) -> bool {
    let u5 = obx5.to_ascii_uppercase();
    let u8 = obx8_flag.to_ascii_uppercase();
    if u5.contains("NOT DETECTED")
        || u5.contains("NEGATIVE")
        || u5.contains("260415000")
        || u5.contains("260385009")
    {
        return false;
    }
    u5.contains("DETECTED")
        || u5.contains("POSITIVE")
        || u5.contains("REACTIVE")
        || u5.contains("260373001")
        || u5.contains("10828004")
        || matches!(u8.as_str(), "A" | "H" | "HH" | "POS" | "DET")
}

/// Parse a batch of pipe-delimited HL7 v2.5.1 `ORU^R01` (ELR) and `ADT` (NSSP syndromic) messages,
/// deduplicating repeat patient specimens via 512-bit `HMAC-SHA256` bigram Bloom filters when `deduplicate_pprl` is true.
pub fn parse_hl7_v2_batch_str(
    raw_hl7: &str,
    pprl_cfg: &PprlConfig,
    deduplicate_pprl: bool,
    source_label: &str,
) -> Result<Hl7v2IngestOutcome, EpiError> {
    // Split raw stream into individual messages starting at MSH
    let mut messages: Vec<Vec<&str>> = Vec::new();
    let mut cur_msg: Vec<&str> = Vec::new();

    for line in raw_hl7.split(['\r', '\n']) {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with("FHS")
            || trimmed.starts_with("BHS")
            || trimmed.starts_with("BTS")
            || trimmed.starts_with("FTS")
        {
            continue;
        }
        if trimmed.starts_with("MSH") {
            if !cur_msg.is_empty() {
                messages.push(cur_msg);
                cur_msg = Vec::new();
            }
        }
        if !cur_msg.is_empty() || trimmed.starts_with("MSH") {
            cur_msg.push(trimmed);
        }
    }
    if !cur_msg.is_empty() {
        messages.push(cur_msg);
    }

    if messages.is_empty() {
        return Err(EpiError::Schema(
            "No HL7 v2.5.1 MSH message segments found in input".to_string(),
        ));
    }

    struct RawParsedMsg {
        ctrl_id: String,
        msg_type: String,
        facility: String,
        iso_date: String,
        abs_day: f64,
        code_str: String,
        detail_str: String,
        is_hit: bool,
        syndrome: NsspSyndromeCategory,
        state_abbr: String,
        demo: PprlDemographicInput,
    }

    let mut raw_parsed = Vec::with_capacity(messages.len());

    for (m_idx, segs) in messages.iter().enumerate() {
        let msh = segs[0];
        let field_sep = msh.chars().nth(3).unwrap_or('|');
        let comp_sep = msh.chars().nth(4).unwrap_or('^');

        let msh_fields: Vec<&str> = msh.split(field_sep).collect();
        // Note: In MSH, MSH-1 is field_sep itself, so msh_fields[1] is MSH-2, msh_fields[3] is MSH-4 (Sending Facility),
        // msh_fields[6] is MSH-7 (Date/Time), msh_fields[8] is MSH-9 (Message Type), msh_fields[9] is MSH-10 (Control ID).
        let facility = msh_fields.get(3).copied().unwrap_or("Facility").to_string();
        let msh_ts = msh_fields.get(6).copied().unwrap_or("20260901");
        let msg_type_raw = msh_fields.get(8).copied().unwrap_or("ORU^R01");
        let msg_type_parts: Vec<&str> = msg_type_raw.split(comp_sep).collect();
        let msg_type = if msg_type_parts.len() >= 2 {
            format!("{}^{}", msg_type_parts[0], msg_type_parts[1])
        } else {
            msg_type_raw.to_string()
        };
        let ctrl_id = msh_fields
            .get(9)
            .copied()
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .unwrap_or_else(|| format!("MSG-{:04}", m_idx + 1));

        let mut family = String::new();
        let mut given = String::new();
        let mut dob = String::new();
        let mut sex = String::new();
        let mut zip = String::new();
        let mut state_abbr = "OR".to_string();
        let mut specimen_ts = String::new();
        let mut obx_codes = Vec::new();
        let mut obx_texts = Vec::new();
        let mut any_elr_pos = false;
        let mut chief_complaint = String::new();
        let mut dg1_code = String::new();

        for &seg in segs.iter().skip(1) {
            let f: Vec<&str> = seg.split(field_sep).collect();
            match f[0] {
                "PID" => {
                    // PID-5 Patient Name (Family^Given), PID-7 DOB, PID-8 Sex, PID-11 Address (Street^^City^State^Zip)
                    if let Some(name_f) = f.get(5) {
                        let np: Vec<&str> = name_f.split(comp_sep).collect();
                        family = np.first().copied().unwrap_or("").to_string();
                        given = np.get(1).copied().unwrap_or("").to_string();
                    }
                    dob = f
                        .get(7)
                        .copied()
                        .unwrap_or("")
                        .chars()
                        .take(8)
                        .collect();
                    sex = f.get(8).copied().unwrap_or("U").to_string();
                    if let Some(addr_f) = f.get(11) {
                        let ap: Vec<&str> = addr_f.split(comp_sep).collect();
                        if let Some(st) = ap.get(3).filter(|s| !s.is_empty()) {
                            state_abbr = (*st).to_string();
                        }
                        if let Some(zp) = ap.get(4) {
                            zip = zp.chars().take(5).collect();
                        }
                    }
                }
                "OBR" => {
                    if let Some(ts) = f.get(7).filter(|s| s.len() >= 8) {
                        specimen_ts = (*ts).to_string();
                    }
                }
                "OBX" => {
                    // OBX-3 Identifier (Code^Text^System), OBX-5 Value, OBX-8 Abnormal Flags
                    let obx3 = f.get(3).copied().unwrap_or("");
                    let c3: Vec<&str> = obx3.split(comp_sep).collect();
                    let code = c3.first().copied().unwrap_or("");
                    let disp = c3.get(1).copied().unwrap_or("");

                    let obx5 = f.get(5).copied().unwrap_or("");
                    let c5: Vec<&str> = obx5.split(comp_sep).collect();
                    let val_disp = if c5.len() >= 2 {
                        format!("{} ({})", c5[1], c5[0])
                    } else {
                        obx5.to_string()
                    };
                    let obx8 = f.get(8).copied().unwrap_or("");

                    // LOINC 8661-1 = Chief complaint (NSSP ADT syndromic segment)
                    if code == "8661-1" || disp.to_ascii_uppercase().contains("CHIEF COMPLAINT") {
                        chief_complaint = val_disp;
                    } else {
                        if !code.is_empty() {
                            obx_codes.push(code.to_string());
                        }
                        obx_texts.push(format!("{disp}: {val_disp}"));
                        if is_elr_result_positive(obx5, obx8) {
                            any_elr_pos = true;
                        }
                    }
                }
                "DG1" => {
                    if let Some(dg_f) = f.get(3) {
                        let dp: Vec<&str> = dg_f.split(comp_sep).collect();
                        dg1_code = dp.first().copied().unwrap_or("").to_string();
                        if chief_complaint.is_empty() {
                            chief_complaint = dp.get(1).copied().unwrap_or(&dg1_code).to_string();
                        }
                    }
                }
                _ => {}
            }
        }

        let eff_ts = if !specimen_ts.is_empty() {
            specimen_ts.as_str()
        } else {
            msh_ts
        };
        let (iso_date, abs_day) =
            hl7_ts_to_iso_and_days(eff_ts).unwrap_or_else(|| ("2026-09-01".to_string(), 0.0));

        let is_adt = msg_type.starts_with("ADT");
        let (code_str, detail_str, is_hit) = if is_adt {
            let cd = if !dg1_code.is_empty() {
                dg1_code.clone()
            } else {
                "8661-1".to_string()
            };
            let det = if !chief_complaint.is_empty() {
                chief_complaint.clone()
            } else {
                obx_texts.join("; ")
            };
            let syn = classify_syndrome(&det, &cd);
            (cd, det, syn != NsspSyndromeCategory::UnspecifiedOther)
        } else {
            let cd = if obx_codes.is_empty() {
                "LOINC-ELR".to_string()
            } else {
                obx_codes.join(",")
            };
            let det = obx_texts.join("; ");
            (cd, det, any_elr_pos)
        };

        let syndrome = classify_syndrome(&detail_str, &code_str);
        let demo = PprlDemographicInput {
            record_id: ctrl_id.clone(),
            source_system: facility.clone(),
            encounter_date: iso_date.clone(),
            encounter_day: abs_day,
            family_name: family,
            given_name: given,
            dob_yyyymmdd: dob,
            sex,
            postal_zip: zip,
        };

        raw_parsed.push(RawParsedMsg {
            ctrl_id,
            msg_type,
            facility,
            iso_date,
            abs_day,
            code_str,
            detail_str,
            is_hit,
            syndrome,
            state_abbr,
            demo,
        });
    }

    // Encode positive / syndromic hit encounters via HMAC-SHA256 Bloom filters and deduplicate
    let hit_indices: Vec<usize> = raw_parsed
        .iter()
        .enumerate()
        .filter_map(|(i, r)| if r.is_hit { Some(i) } else { None })
        .collect();

    let encoded_hits: Vec<EncodedPprlRecord> = hit_indices
        .iter()
        .map(|&idx| encode_pprl_record(&raw_parsed[idx].demo, pprl_cfg))
        .collect();

    let pprl_summary = link_and_deduplicate_pprl(&encoded_hits, pprl_cfg);
    let mut hit_meta_by_raw_idx: HashMap<usize, (String, usize, bool)> = HashMap::new();
    for (pos, &raw_i) in hit_indices.iter().enumerate() {
        let hex_pref = encoded_hits[pos].clk_bloom.to_hex()[0..16].to_string();
        let (_, ent_id, retained) = pprl_summary.entity_cluster_by_record[pos];
        hit_meta_by_raw_idx.insert(raw_i, (hex_pref, ent_id, retained));
    }

    let mut week_counts: BTreeMap<MmwrWeek, f64> = BTreeMap::new();
    let mut syn_map: HashMap<&'static str, usize> = HashMap::new();
    let mut events = Vec::with_capacity(raw_parsed.len());
    let mut dedup_cases = 0usize;

    for (idx, r) in raw_parsed.iter().enumerate() {
        let (bloom_pref, ent_id, pprl_keep) = if let Some(m) = hit_meta_by_raw_idx.get(&idx) {
            m.clone()
        } else {
            let enc = encode_pprl_record(&r.demo, pprl_cfg);
            (enc.clk_bloom.to_hex()[0..16].to_string(), 0, false)
        };

        let count_this = r.is_hit && (!deduplicate_pprl || pprl_keep);
        if count_this {
            dedup_cases += 1;
            *syn_map.entry(r.syndrome.label()).or_insert(0) += 1;
            if let Some(wk) = parse_flexible_date_to_mmwr(&r.iso_date) {
                *week_counts.entry(wk).or_insert(0.0) += 1.0;
            }
        }

        let _ = r.abs_day;
        events.push(ParsedHl7v2Event {
            message_control_id: r.ctrl_id.clone(),
            message_type: r.msg_type.clone(),
            sending_facility: r.facility.clone(),
            event_iso_date: r.iso_date.clone(),
            loinc_or_chief_code: r.code_str.clone(),
            result_or_complaint_text: r.detail_str.clone(),
            is_positive_or_case: r.is_hit,
            syndrome: r.syndrome,
            state_abbr: r.state_abbr.clone(),
            bloom_fingerprint_prefix: bloom_pref,
            linked_entity_id: ent_id,
            retained_after_pprl: count_this,
        });
    }

    let points: Vec<TimeSeriesPoint> = week_counts
        .into_iter()
        .map(|(wk, count)| TimeSeriesPoint {
            date_label: wk.compact_label(),
            epiweek: Some(wk),
            count,
            suppressed: false,
        })
        .collect();

    let mut syndrome_counts: Vec<(String, usize)> = syn_map
        .into_iter()
        .map(|(k, v)| (k.to_string(), v))
        .collect();
    syndrome_counts.sort_by(|a, b| b.1.cmp(&a.1));

    let checksum_hex = EpiStore::compute_dataset_checksum(&points);
    let dataset = TimeSeriesDataset {
        vintage_id: format!("hl7v2_elr_{}_{}", now_unix_secs(), checksum_hex),
        title: format!("HL7 v2.5.1 ELR/NSSP ({source_label})"),
        source_kind: SourceKind::LocalHl7Fhir,
        source_url: format!("hl7v2://{source_label}"),
        ingested_at_utc: format!("epoch:{}", now_unix_secs()),
        default_si_preset: "measles".to_string(),
        checksum_hex,
        points,
    };

    Ok(Hl7v2IngestOutcome {
        dataset,
        total_messages_parsed: raw_parsed.len(),
        positive_or_syndromic_hits: hit_indices.len(),
        deduplicated_incident_cases: dedup_cases,
        syndrome_counts,
        pprl_summary,
        events,
    })
}

pub fn import_hl7_v2_file<P: AsRef<Path>>(
    path: P,
    pprl_cfg: &PprlConfig,
    deduplicate_pprl: bool,
) -> Result<Hl7v2IngestOutcome, EpiError> {
    let pref = path.as_ref();
    let txt = std::fs::read_to_string(pref)?;
    let lbl = pref
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("HL7v2_Batch");
    parse_hl7_v2_batch_str(&txt, pprl_cfg, deduplicate_pprl, lbl)
}

/// Built-in 10-message HL7 v2.5.1 `ORU^R01` (ELR) + `ADT^A04` (CDC NSSP Syndromic) batch stream
/// including repeat PCR + IgM specimens and ED presentation on shared patients (demonstrating Bloom PPRL deduplication)
/// plus a negative specimen that is excluded from case counts.
pub fn builtin_hl7_v2_elr_and_nssp_batch() -> String {
    [
        "MSH|^~\\&|OSP_LAB|OREGON_PHL|EPILEDGER|DOH|20260902091500||ORU^R01^ORU_R01|ELR-9001|P|2.5.1",
        "PID|1||PAT-101^^^OSP||HERNANDEZ^ELENA||19890417|F|||120 NW KINGS^^CORVALLIS^OR^97330",
        "OBR|1||ACC-501|48508-6^Measles virus RNA PCR^LN|||20260901142000",
        "OBX|1|CE|48508-6^Measles virus RNA PCR^LN||260373001^Detected^SCT|||POS|||F",
        // Repeat Emergency Department ADT^A04 visit 2 days later for Elena Hernandez (typo ELENNA) -> deduplicated by Bloom PPRL!
        "MSH|^~\\&|EPIC_ADT|GOOD_SAM_ED|NSSP_BIOSENSE|CDC|20260903181000||ADT^A04^ADT_A01|ADT-9002|P|2.5.1",
        "PID|1||MRN-882^^^GS||HERNANDEZ^ELENNA||19890417|F|||120 NW KINGS^^CORVALLIS^OR^97330",
        "OBX|1|TX|8661-1^Chief complaint^LN||High fever and maculopapular rash||||||F",
        "DG1|1||B05.9^Measles without complication^I10",
        // Distinct Patient 2 (Week 36)
        "MSH|^~\\&|OSP_LAB|OREGON_PHL|EPILEDGER|DOH|20260904113000||ORU^R01^ORU_R01|ELR-9003|P|2.5.1",
        "PID|1||PAT-102^^^OSP||KOWALSKI^JONATHAN||19761123|M|||445 SW 3RD^^CORVALLIS^OR^97333",
        "OBR|1||ACC-509|48508-6^Measles virus RNA PCR^LN|||20260904084500",
        "OBX|1|CE|48508-6^Measles virus RNA PCR^LN||260373001^Detected^SCT|||POS|||F",
        // Confirmatory IgM on Jonathan Kowalski (typo JOHNATHAN) 3 days later -> deduplicated by Bloom PPRL!
        "MSH|^~\\&|QUEST_ELR|REF_LAB|EPILEDGER|DOH|20260907100500||ORU^R01^ORU_R01|ELR-9004|P|2.5.1",
        "PID|1||PAT-992^^^QL||KOWALSKI^JOHNATHAN||19761123|M|||445 SW 3RD^^CORVALLIS^OR^97333",
        "OBR|1||ACC-612|21503-8^Measles virus IgM Ab^LN|||20260907080000",
        "OBX|1|CE|21503-8^Measles virus IgM Ab^LN||10828004^Positive^SCT|||A|||F",
        // Negative PCR rule-out specimen -> ignored by ELR positivity classifier
        "MSH|^~\\&|OSP_LAB|OREGON_PHL|EPILEDGER|DOH|20260908150000||ORU^R01^ORU_R01|ELR-9005|P|2.5.1",
        "PID|1||PAT-105^^^OSP||NGUYEN^MINH||19950214|M|||812 MONROE^^CORVALLIS^OR^97330",
        "OBR|1||ACC-640|48508-6^Measles virus RNA PCR^LN|||20260908103000",
        "OBX|1|CE|48508-6^Measles virus RNA PCR^LN||260415000^Not detected^SCT|||N|||F",
        // Distinct Patient 3 (Week 37 — Pertussis ELR)
        "MSH|^~\\&|OSP_LAB|OREGON_PHL|EPILEDGER|DOH|20260910093000||ORU^R01^ORU_R01|ELR-9006|P|2.5.1",
        "PID|1||PAT-201^^^OSP||CHEN^MARCUS||20180209|M|||90 OAKWAY^^EUGENE^OR^97401",
        "OBR|1||ACC-701|23826-1^Bordetella pertussis DNA PCR^LN|||20260909160000",
        "OBX|1|CE|23826-1^Bordetella pertussis DNA PCR^LN||260373001^Detected^SCT|||POS|||F",
        // Distinct Patient 4 (Week 37 — NSSP ILI ED Visit)
        "MSH|^~\\&|EPIC_ADT|SACRED_HEART_ED|NSSP_BIOSENSE|CDC|20260911211500||ADT^A04^ADT_A01|ADT-9007|P|2.5.1",
        "PID|1||MRN-310^^^SH||CHEN^OLIVIA||20150614|F|||90 OAKWAY^^EUGENE^OR^97401",
        "OBX|1|TX|8661-1^Chief complaint^LN||Paroxysmal cough and fever x 4 days||||||F",
        // Distinct Patient 5 (Week 38 — Measles ELR)
        "MSH|^~\\&|OSP_LAB|OREGON_PHL|EPILEDGER|DOH|20260915134500||ORU^R01^ORU_R01|ELR-9008|P|2.5.1",
        "PID|1||PAT-304^^^OSP||LINDQVIST^ASTRID||19640519|F|||332 CIRCLE BLVD^^CORVALLIS^OR^97330",
        "OBR|1||ACC-788|48508-6^Measles virus RNA PCR^LN|||20260915091500",
        "OBX|1|CE|48508-6^Measles virus RNA PCR^LN||260373001^Detected^SCT|||POS|||F",
        // Distinct Patient 6 (Week 38 — Enteric GI Syndromic)
        "MSH|^~\\&|EPIC_ADT|GOOD_SAM_ED|NSSP_BIOSENSE|CDC|20260917082000||ADT^A04^ADT_A01|ADT-9009|P|2.5.1",
        "PID|1||MRN-419^^^GS||PATEL^DEVIKA||20011203|F|||710 STATE ST^^SALEM^OR^97301",
        "OBX|1|TX|8661-1^Chief complaint^LN||Acute watery diarrhea and vomiting||||||F",
        // Distinct Patient 7 (Week 39 — Measles ELR)
        "MSH|^~\\&|OSP_LAB|OREGON_PHL|EPILEDGER|DOH|20260922160000||ORU^R01^ORU_R01|ELR-9010|P|2.5.1",
        "PID|1||PAT-411^^^OSP||OKAFOR^CHIDI||19930830|M|||500 BURNSIDE^^PORTLAND^OR^97201",
        "OBR|1||ACC-812|48508-6^Measles virus RNA PCR^LN|||20260922110000",
        "OBX|1|CE|48508-6^Measles virus RNA PCR^LN||260373001^Detected^SCT|||POS|||F",
    ]
    .join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hl7_v2_elr_and_nssp_pprl_deduplication() {
        let batch = builtin_hl7_v2_elr_and_nssp_batch();
        let cfg = PprlConfig::default();

        let dedup = parse_hl7_v2_batch_str(&batch, &cfg, true, "TestBatch").expect("parse dedup");
        assert_eq!(dedup.total_messages_parsed, 10);
        // 9 positive/syndromic hits (1 negative PCR excluded), collapsing 2 repeat encounters -> 7 unique incident cases
        assert_eq!(dedup.positive_or_syndromic_hits, 9);
        assert_eq!(dedup.deduplicated_incident_cases, 7);

        let raw = parse_hl7_v2_batch_str(&batch, &cfg, false, "TestBatchRaw").expect("parse raw");
        assert_eq!(raw.deduplicated_incident_cases, 9);
    }
}
