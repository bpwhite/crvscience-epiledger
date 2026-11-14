use epi_core::EpiError;
use serde::{Deserialize, Serialize};

/// Pure-Rust binary reader for SAS Transport Format v5 (`.XPT`), the standard archive format
/// published by CDC NCHS for Continuous NHANES and BRFSS public-use microdata files.
///
/// Reference specification: FDA / SAS TS-140 ("The Record Layout of a Data Set in SAS Transport (XPORT) Format").
/// All records are blocked on 80-byte ASCII cards; numeric values use big-endian IBM System/360
/// hexadecimal floating-point (`IBM370` base-16, excess-64 exponent).

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SasVariableMeta {
    pub name: String,
    pub label: String,
    pub is_char: bool,
    pub length: usize,
    pub offset_in_obs: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum SasXportValue {
    Numeric(f64),
    Missing(char),
    Text(String),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SasXportDataset {
    pub member_name: String,
    pub variables: Vec<SasVariableMeta>,
    pub rows: Vec<Vec<SasXportValue>>,
}

const CARD_LEN: usize = 80;
const LIB_HEADER_PREFIX: &[u8] = b"HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!";
const MEM_HEADER_PREFIX: &[u8] = b"HEADER RECORD*******MEMBER  HEADER RECORD!!!!!!!";
const NAMESTR_HEADER_PREFIX: &[u8] = b"HEADER RECORD*******NAMESTR HEADER RECORD!!!!!!!";
const OBS_HEADER_PREFIX: &[u8] = b"HEADER RECORD*******OBS     HEADER RECORD!!!!!!!";

/// Decode an 8-byte (or 2..=8 byte truncated) big-endian IBM System/360 hexadecimal float into IEEE-754 `f64`,
/// returning `SasXportValue::Missing(tag)` for SAS standard missing `.` (`0x2E`) or special missing `.A`..`.Z` / `._`.
pub fn decode_ibm370_float(bytes: &[u8]) -> SasXportValue {
    if bytes.is_empty() || bytes.len() > 8 {
        return SasXportValue::Missing('.');
    }
    let mut buf = [0u8; 8];
    buf[..bytes.len()].copy_from_slice(bytes);

    // True zero: all 8 bytes are 0x00
    if buf.iter().all(|&b| b == 0) {
        return SasXportValue::Numeric(0.0);
    }

    // SAS missing value marker: byte 0 is '.' (0x2E), '_' (0x5F), or 'A'..='Z', and bytes 1..8 are 0x00
    if buf[1..].iter().all(|&b| b == 0) {
        let b0 = buf[0];
        if b0 == b'.' || b0 == b'_' || b0.is_ascii_uppercase() {
            return SasXportValue::Missing(b0 as char);
        }
    }

    let sign = if (buf[0] & 0x80) != 0 { -1.0 } else { 1.0 };
    let exponent_excess64 = i32::from(buf[0] & 0x7F) - 64;

    // 56-bit unsigned mantissa in bytes 1..8 representing fraction in [1/16, 1)
    let mut mantissa_u64: u64 = 0;
    for &b in &buf[1..] {
        mantissa_u64 = (mantissa_u64 << 8) | u64::from(b);
    }
    if mantissa_u64 == 0 {
        return SasXportValue::Numeric(0.0);
    }

    // Value = sign * (mantissa / 2^56) * 16^(exponent_excess64)
    //       = sign * (mantissa as f64) * 2^(4 * exponent_excess64 - 56)
    let power_of_two = 4 * exponent_excess64 - 56;
    let val = sign * (mantissa_u64 as f64) * 2.0_f64.powi(power_of_two);
    SasXportValue::Numeric(val)
}

/// Parse the first dataset member inside a SAS Transport v5 (`.XPT`) byte buffer (`max_rows` caps memory usage).
pub fn parse_sas_xport_bytes(
    data: &[u8],
    max_rows: usize,
) -> Result<SasXportDataset, EpiError> {
    if data.len() < CARD_LEN * 4 || !data.starts_with(LIB_HEADER_PREFIX) {
        return Err(EpiError::Schema(
            "Input is not a valid SAS XPORT v5 binary file (missing LIBRARY HEADER RECORD)".to_string(),
        ));
    }

    let total_cards = data.len() / CARD_LEN;
    let mut card_idx = 0usize;
    let mut member_name = "MEMBER".to_string();
    let mut namestr_size = 140usize;

    // Scan forward for MEMBER HEADER + NAMESTR HEADER
    while card_idx < total_cards {
        let card = &data[card_idx * CARD_LEN..(card_idx + 1) * CARD_LEN];
        if card.starts_with(MEM_HEADER_PREFIX) {
            // Bytes 74..78 hold NAMESTR record length in ASCII ("0140" on standard v5, "0136" on VAX)
            if let Ok(sz_str) = std::str::from_utf8(&card[74..78]) {
                if let Ok(sz) = sz_str.trim().parse::<usize>() {
                    if sz == 136 || sz == 140 {
                        namestr_size = sz;
                    }
                }
            }
            // Next card after DSCRPTR HEADER holds 8-byte member name at offset 8..16
            if card_idx + 2 < total_cards {
                let ds_card = &data[(card_idx + 2) * CARD_LEN..(card_idx + 3) * CARD_LEN];
                if let Ok(nm) = std::str::from_utf8(&ds_card[8..16]) {
                    let trimmed = nm.trim();
                    if !trimmed.is_empty() {
                        member_name = trimmed.to_string();
                    }
                }
            }
        } else if card.starts_with(NAMESTR_HEADER_PREFIX) {
            break;
        }
        card_idx += 1;
    }

    if card_idx >= total_cards {
        return Err(EpiError::Schema(
            "SAS XPORT file contains no NAMESTR variable descriptor header".to_string(),
        ));
    }

    // Extract variable count from bytes 54..58 of NAMESTR HEADER RECORD
    let namestr_hdr = &data[card_idx * CARD_LEN..(card_idx + 1) * CARD_LEN];
    let nvars: usize = std::str::from_utf8(&namestr_hdr[54..58])
        .unwrap_or("0")
        .trim()
        .parse()
        .map_err(|_| EpiError::Schema("Invalid variable count in NAMESTR header".to_string()))?;

    card_idx += 1;
    let namestr_bytes_needed = nvars * namestr_size;
    let namestr_cards = namestr_bytes_needed.div_ceil(CARD_LEN);
    if card_idx + namestr_cards > total_cards {
        return Err(EpiError::Schema(
            "Truncated NAMESTR descriptor block in SAS XPORT file".to_string(),
        ));
    }

    let namestr_slice =
        &data[card_idx * CARD_LEN..(card_idx + namestr_cards) * CARD_LEN];
    let mut variables = Vec::with_capacity(nvars);
    let mut obs_row_len = 0usize;

    for i in 0..nvars {
        let rec = &namestr_slice[i * namestr_size..(i + 1) * namestr_size];
        let ntype = u16::from_be_bytes([rec[0], rec[1]]);
        let length = u16::from_be_bytes([rec[4], rec[5]]) as usize;
        let name = ascii_trim(&rec[8..16]);
        let label = ascii_trim(&rec[16..56]);
        let npos = u32::from_be_bytes([rec[76], rec[77], rec[78], rec[79]]) as usize;

        obs_row_len = obs_row_len.max(npos + length);
        variables.push(SasVariableMeta {
            name,
            label,
            is_char: ntype == 2,
            length,
            offset_in_obs: npos,
        });
    }

    card_idx += namestr_cards;

    // Advance to OBS HEADER RECORD
    while card_idx < total_cards {
        let card = &data[card_idx * CARD_LEN..(card_idx + 1) * CARD_LEN];
        if card.starts_with(OBS_HEADER_PREFIX) {
            card_idx += 1;
            break;
        }
        card_idx += 1;
    }

    let obs_stream = &data[card_idx * CARD_LEN..];
    let mut rows = Vec::new();
    if obs_row_len == 0 {
        return Ok(SasXportDataset {
            member_name,
            variables,
            rows,
        });
    }

    let mut offset = 0usize;
    while offset + obs_row_len <= obs_stream.len() && rows.len() < max_rows {
        let row_bytes = &obs_stream[offset..offset + obs_row_len];
        // Stop if we hit padding blank cards (all ASCII spaces 0x20) or next HEADER RECORD
        if row_bytes.starts_with(b"HEADER RECORD***") || row_bytes.iter().all(|&b| b == b' ') {
            break;
        }

        let mut parsed_row = Vec::with_capacity(variables.len());
        for var in &variables {
            let cell = &row_bytes[var.offset_in_obs..var.offset_in_obs + var.length];
            if var.is_char {
                parsed_row.push(SasXportValue::Text(ascii_trim(cell)));
            } else {
                parsed_row.push(decode_ibm370_float(cell));
            }
        }
        rows.push(parsed_row);
        offset += obs_row_len;
    }

    Ok(SasXportDataset {
        member_name,
        variables,
        rows,
    })
}

fn ascii_trim(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).trim().to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ibm370_float_conversions() {
        // 1.0 in IBM 370 hex float: +1 * (1/16) * 16^1 -> exp byte 0x41, mantissa 0x10_00_00_00_00_00_00
        let one = decode_ibm370_float(&[0x41, 0x10, 0, 0, 0, 0, 0, 0]);
        assert_eq!(one, SasXportValue::Numeric(1.0));

        // 100.0 in IBM 370 hex float: (100 / 256) * 16^2 -> exp 0x42, mantissa 0x64_00...
        let hundred = decode_ibm370_float(&[0x42, 0x64, 0, 0, 0, 0, 0, 0]);
        assert_eq!(hundred, SasXportValue::Numeric(100.0));

        // -2.5 in IBM 370 hex float: -1 * (2.5 / 16) * 16^1 -> sign|0x41 = 0xC1, 2.5*16 = 40 = 0x28
        let neg = decode_ibm370_float(&[0xC1, 0x28, 0, 0, 0, 0, 0, 0]);
        assert_eq!(neg, SasXportValue::Numeric(-2.5));

        // SAS standard missing '.'
        let miss = decode_ibm370_float(&[b'.', 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(miss, SasXportValue::Missing('.'));
    }
}
