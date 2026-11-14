use serde::{Deserialize, Serialize};
use std::fmt;

/// CDC Morbidity and Mortality Weekly Report (MMWR) epidemiological week.
///
/// Per CSTE / CDC convention, an MMWR week runs Sunday through Saturday (7 days).
/// MMWR Week 1 of year Y is the first week of the year that has at least four days
/// in the calendar year (equivalently, the week containing January 4).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct MmwrWeek {
    pub year: i32,
    pub week: u8,
}

impl MmwrWeek {
    pub fn new(year: i32, week: u8) -> Option<Self> {
        if !(1..=53).contains(&week) || !(1800..=2200).contains(&year) {
            return None;
        }
        Some(Self { year, week })
    }

    /// Convert a Gregorian calendar date `(year, month 1..=12, day 1..=31)` to `MmwrWeek`.
    pub fn from_ymd(year: i32, month: u32, day: u32) -> Option<Self> {
        let jdn = ymd_to_days_since_epoch(year, month, day)?;
        let start_this = mmwr_year_start_days(year);
        let start_next = mmwr_year_start_days(year + 1);

        if jdn >= start_next {
            let wk = ((jdn - start_next) / 7 + 1) as u8;
            Some(Self {
                year: year + 1,
                week: wk,
            })
        } else if jdn >= start_this {
            let wk = ((jdn - start_this) / 7 + 1) as u8;
            Some(Self { year, week: wk })
        } else {
            let start_prev = mmwr_year_start_days(year - 1);
            let wk = ((jdn - start_prev) / 7 + 1) as u8;
            Some(Self {
                year: year - 1,
                week: wk,
            })
        }
    }

    /// Return the Sunday start date `(year, month, day)` of this MMWR week as `YYYY-MM-DD`.
    pub fn sunday_iso_date(&self) -> String {
        let start_days = mmwr_year_start_days(self.year) + i64::from(self.week.saturating_sub(1)) * 7;
        let (y, m, d) = days_since_epoch_to_ymd(start_days);
        format!("{y:04}-{m:02}-{d:02}")
    }

    pub fn compact_label(&self) -> String {
        format!("{}-W{:02}", self.year, self.week)
    }
}

impl fmt::Display for MmwrWeek {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}-W{:02}", self.year, self.week)
    }
}

/// Parse common public health date formats (`YYYY-MM-DD`, `YYYYMMDD`, `YYYY-Www`, `YYYY/MM/DD`) into `MmwrWeek`.
pub fn parse_flexible_date_to_mmwr(raw: &str) -> Option<MmwrWeek> {
    let s = raw.trim();
    if s.is_empty() {
        return None;
    }

    // YYYY-Www or YYYY-W1
    if let Some((y_part, w_part)) = s.split_once("-W").or_else(|| s.split_once("-w")) {
        let yr: i32 = y_part.trim().parse().ok()?;
        let wk: u8 = w_part.trim().parse().ok()?;
        return MmwrWeek::new(yr, wk);
    }

    // YYYYMMDD (e.g. openFDA receivedate)
    if s.len() == 8 && s.bytes().all(|b| b.is_ascii_digit()) {
        let yr: i32 = s[0..4].parse().ok()?;
        let mo: u32 = s[4..6].parse().ok()?;
        let dy: u32 = s[6..8].parse().ok()?;
        return MmwrWeek::from_ymd(yr, mo, dy);
    }

    // YYYY-MM-DD or YYYY/MM/DD (with optional ISO timestamp suffix)
    let date_prefix = s.split(['T', ' ']).next().unwrap_or(s);
    let parts: Vec<&str> = date_prefix.split(['-', '/']).collect();
    if parts.len() == 3 {
        let yr: i32 = parts[0].parse().ok()?;
        let mo: u32 = parts[1].parse().ok()?;
        let dy: u32 = parts[2].parse().ok()?;
        return MmwrWeek::from_ymd(yr, mo, dy);
    }

    None
}

/// Howard Hinnant's civil-from-days / days-from-civil algorithm (days since 1970-01-01 UTC).
fn ymd_to_days_since_epoch(y: i32, m: u32, d: u32) -> Option<i64> {
    if !(1..=12).contains(&m) || !(1..=31).contains(&d) {
        return None;
    }
    let yr = i64::from(y) - if m <= 2 { 1 } else { 0 };
    let era = if yr >= 0 { yr } else { yr - 399 } / 400;
    let yoe = yr - era * 400;
    let mp = i64::from(m) + if m > 2 { -3 } else { 9 };
    let doy = (153 * mp + 2) / 5 + i64::from(d) - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    Some(era * 146097 + doe - 719468)
}

fn days_since_epoch_to_ymd(z: i64) -> (i32, u32, u32) {
    let z_shifted = z + 719468;
    let era = if z_shifted >= 0 {
        z_shifted
    } else {
        z_shifted - 146096
    } / 146097;
    let doe = z_shifted - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let final_y = y + if m <= 2 { 1 } else { 0 };
    (final_y as i32, m as u32, d as u32)
}

/// Day-of-week: Sunday = 0, Monday = 1, ..., Saturday = 6 (1970-01-01 was Thursday = 4).
fn weekday_sunday0(days_since_epoch: i64) -> i64 {
    (days_since_epoch + 4).rem_euclid(7)
}

/// Days since 1970-01-01 of the Sunday starting MMWR Week 1 for `year`
/// (i.e. the Sunday of the week containing January 4 of `year`).
fn mmwr_year_start_days(year: i32) -> i64 {
    let jan4 = ymd_to_days_since_epoch(year, 1, 4).expect("valid Jan 4");
    let dow = weekday_sunday0(jan4);
    jan4 - dow
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mmwr_week_boundaries() {
        // Jan 4, 2026 is a Sunday -> 2026-W01
        let w = MmwrWeek::from_ymd(2026, 1, 4).unwrap();
        assert_eq!(w, MmwrWeek { year: 2026, week: 1 });
        assert_eq!(w.sunday_iso_date(), "2026-01-04");

        // Nov 1, 2026 is a Sunday -> 2026-W44
        let w_nov = parse_flexible_date_to_mmwr("2026-11-01").unwrap();
        assert_eq!(w_nov.year, 2026);
        assert_eq!(w_nov.week, 44);

        // openFDA YYYYMMDD parsing
        let w_fda = parse_flexible_date_to_mmwr("20260825").unwrap();
        assert_eq!(w_fda.year, 2026);
    }
}
