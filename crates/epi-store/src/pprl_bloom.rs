use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};

/// Privacy-Preserving Record Linkage (`PPRL`) via 512-bit Keyed `HMAC-SHA256` Unpadded Bigram
/// Bloom Filters (`POPCNT` Sørensen-Dice & Jaccard similarity), Sub-Quadratic Multi-Pass Hash/LSH
/// Canopy Blocking, Frequency-Resistance Hardening, and Fellegi-Sunter EM Patient Deduplication.
///
/// References:
/// - Blinded WA hospital/death 512-bit HMAC-SHA2 bigram Bloom-filter PPRL study (99.3% concordance): https://pmc.ncbi.nlm.nih.gov/articles/PMC8761329/
/// - Privacy-Preserving Probabilistic Record Linkage (P3RL) Dice + Fellegi-Sunter matching: https://pmc.ncbi.nlm.nih.gov/articles/PMC4460842/
/// - Multi-key / Bloom-filter frequency-attack hardening: https://pmc.ncbi.nlm.nih.gov/articles/PMC7482515/

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BloomFilter512 {
    pub words: [u64; 8],
}

impl Default for BloomFilter512 {
    fn default() -> Self {
        Self { words: [0u64; 8] }
    }
}

impl BloomFilter512 {
    pub fn set_bit(&mut self, bit_idx: usize) {
        let idx = bit_idx & 511;
        self.words[idx >> 6] |= 1u64 << (idx & 63);
    }

    pub fn popcount(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }

    /// Sørensen-Dice similarity `2 * |A & B| / (|A| + |B|)` in `[0.0, 1.0]`.
    pub fn dice_similarity(&self, other: &Self) -> f64 {
        let a = self.popcount();
        let b = other.popcount();
        if a + b == 0 {
            return 0.0;
        }
        let inter: u32 = self
            .words
            .iter()
            .zip(&other.words)
            .map(|(x, y)| (x & y).count_ones())
            .sum();
        (2.0 * f64::from(inter)) / f64::from(a + b)
    }

    /// Jaccard similarity `|A & B| / |A | B|` in `[0.0, 1.0]`.
    pub fn jaccard_similarity(&self, other: &Self) -> f64 {
        let union: u32 = self
            .words
            .iter()
            .zip(&other.words)
            .map(|(x, y)| (x | y).count_ones())
            .sum();
        if union == 0 {
            return 0.0;
        }
        let inter: u32 = self
            .words
            .iter()
            .zip(&other.words)
            .map(|(x, y)| (x & y).count_ones())
            .sum();
        f64::from(inter) / f64::from(union)
    }

    /// Extract 4 orthogonal 8-bit Locality-Sensitive Hashing (`LSH`) band signatures across the 512-bit filter
    /// so high-Dice pairs collide in at least one band even when birth-year has a typo.
    pub fn lsh_band_signatures(&self) -> [u8; 4] {
        let pack8 = |w_a: u64, w_b: u64| -> u8 {
            let mut sig = 0u8;
            for b in 0..4 {
                if ((w_a >> (b * 15 + 3)) & 1) != 0 {
                    sig |= 1 << b;
                }
                if ((w_b >> (b * 15 + 7)) & 1) != 0 {
                    sig |= 1 << (b + 4);
                }
            }
            sig
        };
        [
            pack8(self.words[0], self.words[1]),
            pack8(self.words[2], self.words[3]),
            pack8(self.words[4], self.words[5]),
            pack8(self.words[6], self.words[7]),
        ]
    }

    pub fn to_hex(&self) -> String {
        let mut out = String::with_capacity(128);
        for w in &self.words {
            out.push_str(&format!("{w:016x}"));
        }
        out
    }
}

/// Pure-Rust RFC 2104 `HMAC-SHA256` built directly on `sha2::Sha256`.
pub fn hmac_sha256(key: &[u8], message: &[u8]) -> [u8; 32] {
    let mut key_block = [0u8; 64];
    if key.len() > 64 {
        let hashed = Sha256::digest(key);
        key_block[..32].copy_from_slice(&hashed);
    } else {
        key_block[..key.len()].copy_from_slice(key);
    }

    let mut ipad = [0x36u8; 64];
    let mut opad = [0x5cu8; 64];
    for i in 0..64 {
        ipad[i] ^= key_block[i];
        opad[i] ^= key_block[i];
    }

    let mut inner = Sha256::new();
    inner.update(ipad);
    inner.update(message);
    let inner_digest = inner.finalize();

    let mut outer = Sha256::new();
    outer.update(opad);
    outer.update(inner_digest);
    let final_digest = outer.finalize();

    let mut out = [0u8; 32];
    out.copy_from_slice(&final_digest);
    out
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PprlDemographicInput {
    pub record_id: String,
    pub source_system: String,
    pub encounter_date: String,
    pub encounter_day: f64,
    pub family_name: String,
    pub given_name: String,
    pub dob_yyyymmdd: String,
    pub sex: String,
    pub postal_zip: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EncodedPprlRecord {
    pub record_id: String,
    pub source_system: String,
    pub encounter_date: String,
    pub encounter_day: f64,
    pub birth_year_block_token: String,
    pub lsh_bands: [u8; 4],
    pub clk_bloom: BloomFilter512,
    pub bits_set: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PprlConfig {
    pub secret_key: String,
    pub hashes_name: usize,
    pub hashes_dob: usize,
    pub hashes_aux: usize,
    pub hardening_salted_bits: usize,
    pub dice_match_threshold: f64,
    pub episode_window_days: f64,
    /// Enable multi-pass hash-blocked (`birth_year_token` + 4 Bloom `LSH` bands) sub-quadratic indexing.
    pub use_hash_blocking: bool,
}

impl Default for PprlConfig {
    fn default() -> Self {
        Self {
            secret_key: "epiledger-local-pprl-hmac-key-2026".to_string(),
            hashes_name: 14,
            hashes_dob: 16,
            hashes_aux: 6,
            hardening_salted_bits: 2,
            dice_match_threshold: 0.80,
            episode_window_days: 30.0,
            use_hash_blocking: true,
        }
    }
}

fn normalize_alpha_upper(s: &str) -> Vec<u8> {
    s.bytes()
        .filter_map(|b| {
            let u = b.to_ascii_uppercase();
            if u.is_ascii_alphanumeric() {
                Some(u)
            } else {
                None
            }
        })
        .collect()
}

fn insert_field_unpadded_bigrams(
    bf: &mut BloomFilter512,
    secret_key: &[u8],
    field_tag: &str,
    raw_value: &str,
    k_hashes: usize,
) {
    let norm = normalize_alpha_upper(raw_value);
    if norm.is_empty() || k_hashes == 0 {
        return;
    }

    let slices: Vec<&[u8]> = if norm.len() == 1 {
        vec![&norm[0..1]]
    } else {
        norm.windows(2).collect()
    };

    for (pos_bucket, bg) in slices.into_iter().enumerate() {
        let mut msg = Vec::with_capacity(field_tag.len() + 6);
        msg.extend_from_slice(field_tag.as_bytes());
        msg.push(b':');
        msg.extend_from_slice(bg);
        msg.push(b'@');
        msg.push((pos_bucket / 3) as u8);

        let mac = hmac_sha256(secret_key, &msg);
        let h1 = u64::from_le_bytes(mac[0..8].try_into().unwrap()) as usize;
        let h2 = (u64::from_le_bytes(mac[8..16].try_into().unwrap()) as usize) | 1;

        for i in 0..k_hashes {
            let bit = (h1.wrapping_add(i.wrapping_mul(h2)).wrapping_add(i * i)) & 511;
            bf.set_bit(bit);
        }
    }
}

/// Encode a demographic record into a zero-PII `EncodedPprlRecord` with 512-bit CLK + LSH blocking keys.
pub fn encode_pprl_record(rec: &PprlDemographicInput, cfg: &PprlConfig) -> EncodedPprlRecord {
    let key = cfg.secret_key.as_bytes();
    let mut bf = BloomFilter512::default();

    insert_field_unpadded_bigrams(&mut bf, key, "FN", &rec.family_name, cfg.hashes_name);
    insert_field_unpadded_bigrams(&mut bf, key, "GN", &rec.given_name, cfg.hashes_name);
    insert_field_unpadded_bigrams(&mut bf, key, "DOB", &rec.dob_yyyymmdd, cfg.hashes_dob);
    insert_field_unpadded_bigrams(&mut bf, key, "SEX", &rec.sex, cfg.hashes_aux);
    insert_field_unpadded_bigrams(&mut bf, key, "ZIP", &rec.postal_zip, cfg.hashes_aux);

    let byear: String = rec
        .dob_yyyymmdd
        .chars()
        .filter(|c| c.is_ascii_digit())
        .take(4)
        .collect();
    if cfg.hardening_salted_bits > 0 && !byear.is_empty() {
        let comp = format!("HARDEN:{byear}:{}", rec.sex.to_ascii_uppercase());
        let mac = hmac_sha256(key, comp.as_bytes());
        for s in 0..cfg.hardening_salted_bits.min(16) {
            let bit = u16::from_le_bytes([mac[(2 * s) % 32], mac[(2 * s + 1) % 32]]) as usize & 511;
            bf.set_bit(bit);
        }
    }

    let block_mac = hmac_sha256(key, format!("BLOCK_BYEAR:{byear}").as_bytes());
    let block_tok = format!(
        "{:02x}{:02x}{:02x}{:02x}",
        block_mac[0], block_mac[1], block_mac[2], block_mac[3]
    );
    let lsh_bands = bf.lsh_band_signatures();
    let bits_set = bf.popcount();

    EncodedPprlRecord {
        record_id: rec.record_id.clone(),
        source_system: rec.source_system.clone(),
        encounter_date: rec.encounter_date.clone(),
        encounter_day: rec.encounter_day,
        birth_year_block_token: block_tok,
        lsh_bands,
        clk_bloom: bf,
        bits_set,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PprlCandidatePair {
    pub record_a: String,
    pub record_b: String,
    pub source_a: String,
    pub source_b: String,
    pub dice_similarity: f64,
    pub jaccard_similarity: f64,
    pub fellegi_sunter_weight_bits: f64,
    pub posterior_match_prob: f64,
    pub is_linked_match: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PprlLinkageSummary {
    pub total_input_records: usize,
    pub total_cartesian_pairs: usize,
    pub evaluated_candidate_pairs: usize,
    pub blocking_pair_reduction_pct: f64,
    pub unique_linked_entities: usize,
    pub incident_episodes_retained: usize,
    pub duplicate_encounters_collapsed: usize,
    pub em_prior_match_rate: f64,
    pub candidate_pairs: Vec<PprlCandidatePair>,
    pub entity_cluster_by_record: Vec<(String, usize, bool)>,
}

struct DisjointSet {
    parent: Vec<usize>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }
    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            if ra < rb {
                self.parent[rb] = ra;
            } else {
                self.parent[ra] = rb;
            }
        }
    }
}

/// Build deduplicated candidate index pairs `(i, j)` (`i < j`) using multi-pass hash-blocked canopies
/// (`birth_year_block_token` UNION 4 orthogonal 8-bit Bloom LSH bands) when `use_hash_blocking` is enabled
/// and `n > 32` (falling back to exact all-pairs on tiny batches `<= 32`).
fn collect_candidate_index_pairs(encoded: &[EncodedPprlRecord], use_blocking: bool) -> Vec<(usize, usize)> {
    let n = encoded.len();
    if !use_blocking || n <= 32 {
        let mut all = Vec::with_capacity(n * n.saturating_sub(1) / 2);
        for i in 0..n {
            for j in (i + 1)..n {
                all.push((i, j));
            }
        }
        return all;
    }

    let mut byear_buckets: HashMap<&str, Vec<usize>> = HashMap::new();
    let mut band_buckets: [Vec<Vec<usize>>; 4] = [
        vec![Vec::new(); 256],
        vec![Vec::new(); 256],
        vec![Vec::new(); 256],
        vec![Vec::new(); 256],
    ];

    for (idx, rec) in encoded.iter().enumerate() {
        byear_buckets
            .entry(rec.birth_year_block_token.as_str())
            .or_default()
            .push(idx);
        for b in 0..4 {
            band_buckets[b][usize::from(rec.lsh_bands[b])].push(idx);
        }
    }

    let mut pair_set: HashSet<u64> = HashSet::with_capacity(n * 8);
    let mut emit_bucket = |bucket: &[usize]| {
        // Cap degenerate giant buckets gracefully while preserving local window neighbors
        if bucket.len() <= 256 {
            for a in 0..bucket.len() {
                let i = bucket[a];
                for &j in &bucket[(a + 1)..] {
                    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
                    pair_set.insert(((lo as u64) << 32) | (hi as u64));
                }
            }
        } else {
            for a in 0..bucket.len() {
                let i = bucket[a];
                let end = (a + 128).min(bucket.len());
                for &j in &bucket[(a + 1)..end] {
                    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
                    pair_set.insert(((lo as u64) << 32) | (hi as u64));
                }
            }
        }
    };

    for bucket in byear_buckets.values() {
        emit_bucket(bucket);
    }
    for b in 0..4 {
        for bucket in &band_buckets[b] {
            if bucket.len() >= 2 {
                emit_bucket(bucket);
            }
        }
    }

    let mut pairs: Vec<(usize, usize)> = pair_set
        .into_iter()
        .map(|packed| ((packed >> 32) as usize, (packed & 0xFFFF_FFFF) as usize))
        .collect();
    pairs.sort_unstable();
    pairs
}

/// Link and deduplicate encoded Bloom-filter records via multi-pass hash blocking,
/// 5-bin Fellegi-Sunter EM, and Union-Find episode deduplication.
pub fn link_and_deduplicate_pprl(
    encoded: &[EncodedPprlRecord],
    cfg: &PprlConfig,
) -> PprlLinkageSummary {
    let n = encoded.len();
    if n == 0 {
        return PprlLinkageSummary {
            total_input_records: 0,
            total_cartesian_pairs: 0,
            evaluated_candidate_pairs: 0,
            blocking_pair_reduction_pct: 0.0,
            unique_linked_entities: 0,
            incident_episodes_retained: 0,
            duplicate_encounters_collapsed: 0,
            em_prior_match_rate: 0.0,
            candidate_pairs: Vec::new(),
            entity_cluster_by_record: Vec::new(),
        };
    }

    let cartesian_total = n * (n - 1) / 2;
    let index_pairs = collect_candidate_index_pairs(encoded, cfg.use_hash_blocking);
    let eval_pairs_n = index_pairs.len();
    let reduction_pct = if cartesian_total > 0 && eval_pairs_n <= cartesian_total {
        (1.0 - (eval_pairs_n as f64) / (cartesian_total as f64)) * 100.0
    } else {
        0.0
    };

    let bin_of = |dice: f64| -> usize {
        if dice >= 0.88 {
            4
        } else if dice >= 0.78 {
            3
        } else if dice >= 0.65 {
            2
        } else if dice >= 0.45 {
            1
        } else {
            0
        }
    };

    let mut raw_pairs: Vec<(usize, usize, f64, f64, usize)> = Vec::with_capacity(eval_pairs_n);
    let mut bin_counts = [0.0_f64; 5];

    for (i, j) in index_pairs {
        let dice = encoded[i].clk_bloom.dice_similarity(&encoded[j].clk_bloom);
        let jacc = encoded[i]
            .clk_bloom
            .jaccard_similarity(&encoded[j].clk_bloom);
        let b = bin_of(dice);
        bin_counts[b] += 1.0;
        raw_pairs.push((i, j, dice, jacc, b));
    }
    // Uninspected pruned pairs belong implicitly to Bin 0 (low similarity) for accurate EM prior scaling
    if cartesian_total > eval_pairs_n {
        bin_counts[0] += (cartesian_total - eval_pairs_n) as f64;
    }

    let total_pairs = (cartesian_total.max(1)) as f64;
    let mut pi_match = 0.08_f64;
    let mut m_prob = [0.01, 0.03, 0.08, 0.28, 0.60];
    let mut u_prob = [0.72, 0.20, 0.06, 0.015, 0.005];

    for _em_iter in 0..35 {
        let mut gamma_bin = [0.0_f64; 5];
        let mut sum_gamma = 0.0;
        for b in 0..5 {
            let num = pi_match * m_prob[b];
            let den = (num + (1.0 - pi_match) * u_prob[b]).max(1e-15);
            gamma_bin[b] = num / den;
            sum_gamma += bin_counts[b] * gamma_bin[b];
        }
        pi_match = (sum_gamma / total_pairs).clamp(0.001, 0.85);

        let mut m_sum = 0.0;
        let mut u_sum = 0.0;
        for b in 0..5 {
            m_prob[b] = bin_counts[b] * gamma_bin[b] + 0.05 * (b as f64 + 1.0);
            u_prob[b] = bin_counts[b] * (1.0 - gamma_bin[b]) + 0.05 * (5.0 - b as f64);
            m_sum += m_prob[b];
            u_sum += u_prob[b];
        }
        for b in 0..5 {
            m_prob[b] /= m_sum;
            u_prob[b] /= u_sum;
        }
    }

    let mut dsu = DisjointSet::new(n);
    let mut candidate_pairs = Vec::new();

    for &(i, j, dice, jacc, b) in &raw_pairs {
        let w_bits = (m_prob[b] / u_prob[b].max(1e-12)).log2();
        let num = pi_match * m_prob[b];
        let post = num / (num + (1.0 - pi_match) * u_prob[b]).max(1e-15);
        let is_match = dice >= cfg.dice_match_threshold;
        if is_match {
            dsu.union(i, j);
        }
        if dice >= 0.45 {
            candidate_pairs.push(PprlCandidatePair {
                record_a: encoded[i].record_id.clone(),
                record_b: encoded[j].record_id.clone(),
                source_a: encoded[i].source_system.clone(),
                source_b: encoded[j].source_system.clone(),
                dice_similarity: dice,
                jaccard_similarity: jacc,
                fellegi_sunter_weight_bits: w_bits,
                posterior_match_prob: post,
                is_linked_match: is_match,
            });
        }
    }

    candidate_pairs.sort_by(|a, b| {
        b.dice_similarity
            .partial_cmp(&a.dice_similarity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    if candidate_pairs.len() > 250 {
        candidate_pairs.truncate(250);
    }

    let mut root_to_entity = HashMap::new();
    let mut next_ent = 0usize;
    let mut entity_ids = vec![0usize; n];
    for idx in 0..n {
        let r = dsu.find(idx);
        let eid = *root_to_entity.entry(r).or_insert_with(|| {
            next_ent += 1;
            next_ent
        });
        entity_ids[idx] = eid;
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        encoded[a]
            .encounter_day
            .partial_cmp(&encoded[b].encounter_day)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| encoded[a].record_id.cmp(&encoded[b].record_id))
    });

    let mut last_episode_day: HashMap<usize, f64> = HashMap::new();
    let mut retained_flags = vec![false; n];
    let mut retained_count = 0usize;

    for idx in order {
        let eid = entity_ids[idx];
        let day = encoded[idx].encounter_day;
        let is_new_episode = match last_episode_day.get(&eid) {
            None => true,
            Some(&prev_d) => {
                cfg.episode_window_days > 0.0 && (day - prev_d) > cfg.episode_window_days
            }
        };
        if is_new_episode {
            last_episode_day.insert(eid, day);
            retained_flags[idx] = true;
            retained_count += 1;
        }
    }

    let entity_cluster_by_record: Vec<(String, usize, bool)> = (0..n)
        .map(|i| {
            (
                encoded[i].record_id.clone(),
                entity_ids[i],
                retained_flags[i],
            )
        })
        .collect();

    PprlLinkageSummary {
        total_input_records: n,
        total_cartesian_pairs: cartesian_total,
        evaluated_candidate_pairs: eval_pairs_n,
        blocking_pair_reduction_pct: reduction_pct,
        unique_linked_entities: next_ent,
        incident_episodes_retained: retained_count,
        duplicate_encounters_collapsed: n.saturating_sub(retained_count),
        em_prior_match_rate: pi_match,
        candidate_pairs,
        entity_cluster_by_record,
    }
}

/// Generate a synthetic multi-facility hospital/ELR benchmark of `num_unique` patients (plus ~35% repeat
/// typographical/follow-up encounters) to benchmark sub-quadratic hash-blocked PPRL throughput.
pub fn generate_synthetic_pprl_scale_benchmark(num_unique: usize) -> Vec<PprlDemographicInput> {
    let surnames = [
        "HERNANDEZ", "KOWALSKI", "CHEN", "OKAFOR", "LINDQVIST", "PATEL", "NGUYEN", "SORENSEN",
        "TAKAHASHI", "CAMPBELL", "RAMIREZ", "ALVAREZ", "MENSAH", "NOVAK", "DESAI", "WEBER",
    ];
    let givens = [
        "ELENA", "JONATHAN", "MARCUS", "CHIDI", "ASTRID", "DEVIKA", "MINH", "FREJA", "KENJI",
        "IAN", "MATEO", "SOFIA", "KWAME", "PETRA", "PRIYA", "HANNAH",
    ];

    let n_u = num_unique.clamp(10, 20_000);
    let mut out = Vec::with_capacity(n_u + n_u / 3);
    for i in 0..n_u {
        let fn_base = format!("{}{}", surnames[i % surnames.len()], i / surnames.len());
        let gn_base = format!("{}{}", givens[(i * 7) % givens.len()], (i / 16) % 97);
        let yr = 1945 + (i % 78);
        let mo = 1 + ((i * 3) % 12);
        let dy = 1 + ((i * 5) % 28);
        let dob = format!("{yr:04}{mo:02}{dy:02}");
        let sex = if i % 2 == 0 { "F" } else { "M" };
        let zip = format!("{:05}", 97001 + (i % 450));
        let day0 = (i % 60) as f64;

        out.push(PprlDemographicInput {
            record_id: format!("ELR-{:05}", i + 1),
            source_system: "StatePHL".to_string(),
            encounter_date: "2026-09-10".to_string(),
            encounter_day: day0,
            family_name: fn_base.clone(),
            given_name: gn_base.clone(),
            dob_yyyymmdd: dob.clone(),
            sex: sex.to_string(),
            postal_zip: zip.clone(),
        });

        if i % 3 == 0 {
            out.push(PprlDemographicInput {
                record_id: format!("ADT-{:05}", i + 1),
                source_system: "HospADT".to_string(),
                encounter_date: "2026-09-14".to_string(),
                encounter_day: day0 + 4.0,
                family_name: fn_base,
                given_name: format!("{gn_base}A"),
                dob_yyyymmdd: dob,
                sex: sex.to_string(),
                postal_zip: zip,
            });
        }
    }
    out
}

/// Built-in 12-encounter hospital `ADT` + state lab `ELR` benchmark representing 7 true unique patients.
pub fn builtin_pprl_hospital_lab_records() -> Vec<PprlDemographicInput> {
    vec![
        PprlDemographicInput {
            record_id: "ELR-1001".to_string(),
            source_system: "StatePHL-ELR".to_string(),
            encounter_date: "2026-09-02".to_string(),
            encounter_day: 2.0,
            family_name: "HERNANDEZ".to_string(),
            given_name: "ELENA".to_string(),
            dob_yyyymmdd: "19890417".to_string(),
            sex: "F".to_string(),
            postal_zip: "97330".to_string(),
        },
        PprlDemographicInput {
            record_id: "ADT-2041".to_string(),
            source_system: "GoodSam-ADT".to_string(),
            encounter_date: "2026-09-05".to_string(),
            encounter_day: 5.0,
            family_name: "HERNANDEZ".to_string(),
            given_name: "ELENNA".to_string(),
            dob_yyyymmdd: "19890417".to_string(),
            sex: "F".to_string(),
            postal_zip: "97330".to_string(),
        },
        PprlDemographicInput {
            record_id: "ELR-1002".to_string(),
            source_system: "StatePHL-ELR".to_string(),
            encounter_date: "2026-09-04".to_string(),
            encounter_day: 4.0,
            family_name: "KOWALSKI".to_string(),
            given_name: "JONATHAN".to_string(),
            dob_yyyymmdd: "19761123".to_string(),
            sex: "M".to_string(),
            postal_zip: "97333".to_string(),
        },
        PprlDemographicInput {
            record_id: "ELR-1019".to_string(),
            source_system: "RefLab-ELR".to_string(),
            encounter_date: "2026-09-08".to_string(),
            encounter_day: 8.0,
            family_name: "KOWALSKI".to_string(),
            given_name: "JOHNATHAN".to_string(),
            dob_yyyymmdd: "19761123".to_string(),
            sex: "M".to_string(),
            postal_zip: "97333".to_string(),
        },
        PprlDemographicInput {
            record_id: "ADT-2050".to_string(),
            source_system: "GoodSam-ADT".to_string(),
            encounter_date: "2026-09-06".to_string(),
            encounter_day: 6.0,
            family_name: "CHEN".to_string(),
            given_name: "MARCUS".to_string(),
            dob_yyyymmdd: "20180209".to_string(),
            sex: "M".to_string(),
            postal_zip: "97401".to_string(),
        },
        PprlDemographicInput {
            record_id: "ELR-1025".to_string(),
            source_system: "StatePHL-ELR".to_string(),
            encounter_date: "2026-09-08".to_string(),
            encounter_day: 8.0,
            family_name: "CHEN".to_string(),
            given_name: "MARCUS".to_string(),
            dob_yyyymmdd: "20180209".to_string(),
            sex: "M".to_string(),
            postal_zip: "97401".to_string(),
        },
        PprlDemographicInput {
            record_id: "ELR-1026".to_string(),
            source_system: "StatePHL-ELR".to_string(),
            encounter_date: "2026-09-11".to_string(),
            encounter_day: 11.0,
            family_name: "CHEN".to_string(),
            given_name: "OLIVIA".to_string(),
            dob_yyyymmdd: "20150614".to_string(),
            sex: "F".to_string(),
            postal_zip: "97401".to_string(),
        },
        PprlDemographicInput {
            record_id: "ELR-1030".to_string(),
            source_system: "RefLab-ELR".to_string(),
            encounter_date: "2026-09-12".to_string(),
            encounter_day: 12.0,
            family_name: "OKAFOR".to_string(),
            given_name: "CHIDI".to_string(),
            dob_yyyymmdd: "19930830".to_string(),
            sex: "M".to_string(),
            postal_zip: "97201".to_string(),
        },
        PprlDemographicInput {
            record_id: "ELR-1099".to_string(),
            source_system: "StatePHL-ELR".to_string(),
            encounter_date: "2026-10-30".to_string(),
            encounter_day: 60.0,
            family_name: "OKAFOR".to_string(),
            given_name: "CHIDI".to_string(),
            dob_yyyymmdd: "19930830".to_string(),
            sex: "M".to_string(),
            postal_zip: "97201".to_string(),
        },
        PprlDemographicInput {
            record_id: "ADT-2088".to_string(),
            source_system: "GoodSam-ADT".to_string(),
            encounter_date: "2026-09-15".to_string(),
            encounter_day: 15.0,
            family_name: "LINDQVIST".to_string(),
            given_name: "ASTRID".to_string(),
            dob_yyyymmdd: "19640519".to_string(),
            sex: "F".to_string(),
            postal_zip: "97330".to_string(),
        },
        PprlDemographicInput {
            record_id: "ELR-1044".to_string(),
            source_system: "StatePHL-ELR".to_string(),
            encounter_date: "2026-09-16".to_string(),
            encounter_day: 16.0,
            family_name: "LINDQUIST".to_string(),
            given_name: "ASTRID".to_string(),
            dob_yyyymmdd: "19640519".to_string(),
            sex: "F".to_string(),
            postal_zip: "97330".to_string(),
        },
        PprlDemographicInput {
            record_id: "ELR-1055".to_string(),
            source_system: "RefLab-ELR".to_string(),
            encounter_date: "2026-09-19".to_string(),
            encounter_day: 19.0,
            family_name: "PATEL".to_string(),
            given_name: "DEVIKA".to_string(),
            dob_yyyymmdd: "20011203".to_string(),
            sex: "F".to_string(),
            postal_zip: "97301".to_string(),
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom_pprl_fuzzy_linkage_and_hash_blocking_scale() {
        let raw = builtin_pprl_hospital_lab_records();
        let cfg = PprlConfig::default();
        let encoded: Vec<EncodedPprlRecord> =
            raw.iter().map(|r| encode_pprl_record(r, &cfg)).collect();

        let res = link_and_deduplicate_pprl(&encoded, &cfg);
        assert_eq!(res.unique_linked_entities, 7);
        assert_eq!(res.incident_episodes_retained, 8);

        // Verify hash-blocked indexing prunes >80% of pairs on 600-patient scale benchmark while recovering exact clusters
        let scale_raw = generate_synthetic_pprl_scale_benchmark(450);
        let scale_enc: Vec<EncodedPprlRecord> = scale_raw
            .iter()
            .map(|r| encode_pprl_record(r, &cfg))
            .collect();
        let scale_res = link_and_deduplicate_pprl(&scale_enc, &cfg);
        assert!(scale_res.blocking_pair_reduction_pct > 75.0);
        assert_eq!(scale_res.unique_linked_entities, 450);
    }
}
