// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

use epi_core::{
    parse_flexible_date_to_mmwr, EpiError, SourceKind, TimeSeriesDataset, TimeSeriesPoint,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Deterministic CDC NSSP ESSENCE-style Concatenated Chief-Complaint & Discharge-Diagnosis (`CCDD`)
/// syndrome classification with clinical context-window negation (`NegEx`) and boolean `OR` / `ANDNOT` queries.
///
/// References:
/// - Evaluation and Revision of NSSP ESSENCE Syndromic Definitions: https://pmc.ncbi.nlm.nih.gov/articles/PMC11556616/
/// - Chief-Complaint and ICD-10-CM Exclusion Validation in NSSP: https://pmc.ncbi.nlm.nih.gov/articles/PMC11426074/
/// - All-Hazards Emergency Department Syndromic Surveillance: https://pmc.ncbi.nlm.nih.gov/articles/PMC5676514/

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum EssenceSyndromeId {
    IliRespiratory,
    CovidLikeIllness,
    AcuteGiEnteric,
    SuspectedOverdose,
    HeatRelatedIllness,
    FeverRashExanthem,
}

impl EssenceSyndromeId {
    pub const ALL: [Self; 6] = [
        Self::IliRespiratory,
        Self::CovidLikeIllness,
        Self::AcuteGiEnteric,
        Self::SuspectedOverdose,
        Self::HeatRelatedIllness,
        Self::FeverRashExanthem,
    ];

    pub fn code(self) -> &'static str {
        match self {
            Self::IliRespiratory => "ILI_NSSP_V2",
            Self::CovidLikeIllness => "CLI_NSSP_V2",
            Self::AcuteGiEnteric => "GI_ENTERIC_V1",
            Self::SuspectedOverdose => "DOSE_OVERDOSE_V2",
            Self::HeatRelatedIllness => "HRI_HEAT_V2",
            Self::FeverRashExanthem => "FEVER_RASH_V1",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::IliRespiratory => "Influenza-Like Illness (ILI)",
            Self::CovidLikeIllness => "COVID-Like Illness (CLI)",
            Self::AcuteGiEnteric => "Acute GI / Enteric",
            Self::SuspectedOverdose => "Suspected Unintentional Overdose (DOSE)",
            Self::HeatRelatedIllness => "Environmental Heat-Related Illness (HRI)",
            Self::FeverRashExanthem => "Fever + Rash / Exanthem",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdTriageEncounter {
    pub visit_id: String,
    pub visit_iso_date: String,
    pub facility_id: String,
    pub patient_age: u8,
    pub chief_complaint_text: String,
    /// Semicolon- or space-separated ICD-10-CM discharge diagnosis codes (may be empty at initial triage).
    pub icd10_cm_codes: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassifiedEdEncounter {
    pub visit_id: String,
    pub visit_iso_date: String,
    pub facility_id: String,
    pub normalized_cc_text: String,
    pub matched_syndromes: Vec<EssenceSyndromeId>,
    pub matched_inclusion_terms: Vec<String>,
    pub suppressed_negated_terms: Vec<String>,
    pub custom_query_hit: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyndromeSummaryRow {
    pub syndrome: EssenceSyndromeId,
    pub code: String,
    pub label: String,
    pub matched_visits: usize,
    pub pct_of_ed_visits: f64,
    pub negated_false_positives_blocked: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EssenceBatchClassification {
    pub total_encounters: usize,
    pub encounters_matching_any_syndrome: usize,
    pub total_negations_blocked: usize,
    pub custom_query_expr: String,
    pub custom_query_matches: usize,
    pub syndrome_summaries: Vec<SyndromeSummaryRow>,
    pub classified_encounters: Vec<ClassifiedEdEncounter>,
}

const PRE_NEGATION_TRIGGERS: &[&str] = &[
    "denies",
    "denied",
    "deny",
    "no",
    "not",
    "without",
    "negative",
    "neg",
    "rule out",
    "r o",
    "ruled out for",
    "ruled out of",
    "free of",
    "no evidence of",
    "absent",
    "resolved",
    "hx of",
    "history of",
];

const POST_NEGATION_TRIGGERS: &[&str] = &[
    "ruled out",
    "rule out",
    "r o",
    "negative",
    "neg",
    "not detected",
    "undetected",
    "no evidence",
    "absent",
    "unlikely",
    "excluded",
    "resolved",
    "not present",
    "not found",
];

const POSITIVE_POLARITY_MARKERS: &[&str] = &["positive", "pos", "present", "confirmed"];

const CLAUSE_BREAKERS: &[&str] = &[
    ";",
    ".",
    "but",
    "however",
    "except",
    "now",
    "presents",
    "presenting",
    "complains",
    "c o",
];

fn is_comma_clause_boundary_after(after: &[&str]) -> bool {
    let is_boundary = |tok: &str| {
        matches!(tok, ";" | "." | ",") || CLAUSE_BREAKERS.contains(&tok)
    };
    let Some(&first) = after.first() else {
        return false;
    };
    if is_boundary(first) {
        return false;
    }

    // Unambiguous pre-negation intros
    if matches!(first, "denies" | "denied" | "deny" | "without") {
        return true;
    }
    if after.starts_with(&["free", "of"])
        || after.starts_with(&["hx", "of"])
        || after.starts_with(&["history", "of"])
        || after.starts_with(&["no", "evidence", "of"])
        || after.starts_with(&["ruled", "out", "for"])
        || after.starts_with(&["ruled", "out", "of"])
    {
        return true;
    }

    // Expanded contraction negation ("doesn't", "don't", "can't", etc.)
    if after.len() >= 2
        && after[1] == "not"
        && matches!(
            first,
            "do" | "does"
                | "did"
                | "can"
                | "is"
                | "are"
                | "was"
                | "were"
                | "has"
                | "have"
                | "will"
                | "would"
                | "could"
                | "should"
        )
    {
        return true;
    }

    // "no" pre-negates unless it is trailing post-negation "no evidence"
    if first == "no" {
        if after.starts_with(&["no", "evidence"]) {
            return after.len() > 2 && !is_boundary(after[2]);
        }
        return true;
    }

    // "not" pre-negates (and isolates "fever, not ruled out") unless it is trailing
    // post-negation "not detected" / "not present" / "not found"
    if first == "not" {
        if after.len() >= 2
            && matches!(after[1], "detected" | "present" | "found")
            && (after.len() == 2 || is_boundary(after[2]))
        {
            return false;
        }
        return true;
    }

    // Dual-use 2-token pre-negations ("rule out X", "r o X") and infinitive "to rule out X" / "to r o X" when followed by content
    if after.starts_with(&["rule", "out"]) || after.starts_with(&["r", "o"]) {
        return after.len() > 2 && !is_boundary(after[2]);
    }
    if after.starts_with(&["to", "rule", "out"]) || after.starts_with(&["to", "r", "o"]) {
        return after.len() > 3 && !is_boundary(after[3]);
    }

    // Dual-use 1-token triggers ("negative", "neg", "absent", "resolved") act as pre-negations
    // when introducing a following target ("negative for fever", "absent cough"), and remain
    // post-negations when terminating a clause ("covid, negative" / "fever, absent").
    if matches!(first, "negative" | "neg") {
        return after.len() >= 2
            && (matches!(after[1], "for" | "of")
                || (!is_boundary(after[1]) && !POSITIVE_POLARITY_MARKERS.contains(&after[1])));
    }
    if matches!(first, "absent" | "resolved") {
        return after.len() >= 2
            && !is_boundary(after[1])
            && !POSITIVE_POLARITY_MARKERS.contains(&after[1]);
    }

    false
}

/// Normalize ED triage shorthand and punctuation while preserving clause boundary markers (`;` and `.`),
/// keeping decimal points between ASCII digits (e.g., `102.4`) intact, preserving commas before
/// pre-negation clauses as `;`, collapsing `c/o` (`c o`) to `;`, and expanding common contractions.
pub fn normalize_chief_complaint(raw: &str) -> String {
    let chars: Vec<char> = raw.chars().collect();
    let mut lowered = String::with_capacity(raw.len() + 16);
    for (idx, &ch) in chars.iter().enumerate() {
        if ch.is_ascii_alphanumeric() || ch.is_whitespace() {
            lowered.push(ch.to_ascii_lowercase());
        } else if ch == '.' {
            let between_digits = idx > 0
                && idx + 1 < chars.len()
                && chars[idx - 1].is_ascii_digit()
                && chars[idx + 1].is_ascii_digit();
            if between_digits {
                lowered.push('.');
            } else {
                lowered.push_str(" . ");
            }
        } else if ch == ';' {
            lowered.push_str(" ; ");
        } else if ch == ',' {
            lowered.push_str(" , ");
        } else if matches!(ch, '/' | '-' | ':' | '(' | ')' | '[' | ']') {
            lowered.push(' ');
        }
    }

    // Expand standard ED triage abbreviations and common English contractions token-by-token
    let mut expanded_tokens: Vec<&str> = Vec::new();
    let raw_toks: Vec<&str> = lowered.split_whitespace().collect();
    let mut i = 0;
    while i < raw_toks.len() {
        if i + 2 < raw_toks.len()
            && raw_toks[i] == "n"
            && raw_toks[i + 1] == "v"
            && raw_toks[i + 2] == "d"
        {
            expanded_tokens.extend(["nausea", "vomiting", "diarrhea"]);
            i += 3;
            continue;
        }
        if i + 1 < raw_toks.len() && raw_toks[i] == "n" && raw_toks[i + 1] == "v" {
            expanded_tokens.extend(["nausea", "vomiting"]);
            i += 2;
            continue;
        }
        if i + 1 < raw_toks.len() && raw_toks[i] == "c" && raw_toks[i + 1] == "o" {
            expanded_tokens.push(";");
            i += 2;
            continue;
        }
        match raw_toks[i] {
            "sob" => expanded_tokens.extend(["shortness", "of", "breath"]),
            "uri" => expanded_tokens.extend(["upper", "respiratory", "infection"]),
            "ili" => expanded_tokens.extend(["influenza", "like", "illness"]),
            "od" => expanded_tokens.push("overdose"),
            "ha" => expanded_tokens.push("headache"),
            "fev" | "tmax" => expanded_tokens.push("fever"),
            "vomit" | "emesis" => expanded_tokens.push("vomiting"),
            "dont" => expanded_tokens.extend(["do", "not"]),
            "doesnt" => expanded_tokens.extend(["does", "not"]),
            "didnt" => expanded_tokens.extend(["did", "not"]),
            "cant" => expanded_tokens.extend(["can", "not"]),
            "isnt" => expanded_tokens.extend(["is", "not"]),
            "arent" => expanded_tokens.extend(["are", "not"]),
            "wasnt" => expanded_tokens.extend(["was", "not"]),
            "werent" => expanded_tokens.extend(["were", "not"]),
            "hasnt" => expanded_tokens.extend(["has", "not"]),
            "havent" => expanded_tokens.extend(["have", "not"]),
            "wont" => expanded_tokens.extend(["will", "not"]),
            "wouldnt" => expanded_tokens.extend(["would", "not"]),
            "couldnt" => expanded_tokens.extend(["could", "not"]),
            "shouldnt" => expanded_tokens.extend(["should", "not"]),
            other => expanded_tokens.push(other),
        }
        i += 1;
    }

    // Resolve `,` tokens: preserve as `;` before pre-negation clauses, drop before plain coordinate items.
    let mut resolved_tokens: Vec<&str> = Vec::with_capacity(expanded_tokens.len());
    for (idx, &tok) in expanded_tokens.iter().enumerate() {
        if tok == "," {
            if is_comma_clause_boundary_after(&expanded_tokens[idx + 1..]) {
                if resolved_tokens.last().copied() != Some(";") {
                    resolved_tokens.push(";");
                }
            }
        } else {
            resolved_tokens.push(tok);
        }
    }
    resolved_tokens.join(" ")
}

fn slice_contains_phrase_toks(tokens: &[&str], phrase: &str) -> bool {
    let trig_toks: Vec<&str> = phrase.split_whitespace().collect();
    if trig_toks.is_empty() || tokens.len() < trig_toks.len() {
        return false;
    }
    tokens
        .windows(trig_toks.len())
        .any(|w| w == trig_toks.as_slice())
}

/// Check if a candidate symptom phrase appears affirmatively (un-negated) or only under negation.
/// Evaluates both pre-negation triggers (preceding tokens up to clause boundary / prior post-negation)
/// and NegEx post-negation triggers (trailing tokens up to clause boundary or positive markers).
/// Returns `(affirmative_hit, negated_hit)`.
pub fn match_affirmative_or_negated(
    normalized_cc: &str,
    phrase: &str,
    neg_window_tokens: usize,
) -> (bool, bool) {
    let target_toks: Vec<&str> = phrase.split_whitespace().collect();
    if target_toks.is_empty() {
        return (false, false);
    }
    let cc_toks: Vec<&str> = normalized_cc.split_whitespace().collect();
    if cc_toks.len() < target_toks.len() {
        return (false, false);
    }

    let win = neg_window_tokens.max(1);
    let mut found_affirmative = false;
    let mut found_negated = false;

    for pos in 0..=(cc_toks.len() - target_toks.len()) {
        if cc_toks[pos..pos + target_toks.len()] == target_toks[..] {
            // 1. Pre-negation lookback up to `win` tokens, stopping at clause breakers
            // or unnegated positive polarity markers (excluding "not present" / "not positive")
            let mut clause_start = 0usize;
            for back in (0..pos).rev() {
                let tok = cc_toks[back];
                let is_positive_breaker = POSITIVE_POLARITY_MARKERS.contains(&tok)
                    && !(back > 0 && matches!(cc_toks[back - 1], "not" | "no" | "never"));
                if CLAUSE_BREAKERS.contains(&tok)
                    || (tok == "o" && back > 0 && cc_toks[back - 1] == "c")
                    || is_positive_breaker
                {
                    clause_start = back + 1;
                    break;
                }
            }
            let win_start = pos.saturating_sub(win);
            let active_slice_start = win_start.max(clause_start);

            // Prevent a prior symptom's post-negation phrase from bleeding forward across
            // punctuation-free clauses into `pos`, while preserving true pre-negation intros
            // (e.g. "no evidence of measles", "negative for fever", clause-initial "rule out fever",
            // or "to rule out covid").
            let mut prior_post_cut = active_slice_start;
            for trig in POST_NEGATION_TRIGGERS {
                let trig_toks: Vec<&str> = trig.split_whitespace().collect();
                let tlen = trig_toks.len();
                if tlen == 0 || pos < clause_start + tlen {
                    continue;
                }
                let is_dual_use = matches!(
                    *trig,
                    "absent"
                        | "resolved"
                        | "rule out"
                        | "r o"
                        | "negative"
                        | "neg"
                        | "no evidence"
                );
                for idx in clause_start..=(pos - tlen) {
                    if cc_toks[idx..idx + tlen] == trig_toks[..] {
                        let trig_end = idx + tlen;
                        if matches!(cc_toks.get(trig_end).copied(), Some("for" | "of")) {
                            continue;
                        }
                        if is_dual_use
                            && (idx == clause_start
                                || (matches!(*trig, "rule out" | "r o")
                                    && idx > clause_start
                                    && cc_toks[idx - 1] == "to"))
                        {
                            continue;
                        }
                        prior_post_cut = prior_post_cut.max(trig_end);
                    }
                }
            }

            let pre_slice_start = active_slice_start.max(prior_post_cut);
            let pre_tokens = &cc_toks[pre_slice_start..pos];
            let is_pre_neg = PRE_NEGATION_TRIGGERS
                .iter()
                .any(|trig| slice_contains_phrase_toks(pre_tokens, trig));

            // 2. Post-negation forward scan up to `win` tokens, stopping at clause breakers,
            // positive markers (except "not present" which is a post-negation trigger),
            // or pre-negation intros ("negative for", "no evidence of", "to rule out", etc.).
            let post_start = pos + target_toks.len();
            let post_max = (post_start + win).min(cc_toks.len());
            let mut post_end = post_max;
            for fwd in post_start..post_max {
                let tok = cc_toks[fwd];
                let next_is_prep = matches!(cc_toks.get(fwd + 1).copied(), Some("for" | "of"));
                let is_to_rule_out = tok == "to"
                    && ((cc_toks.get(fwd + 1).copied() == Some("rule")
                        && cc_toks.get(fwd + 2).copied() == Some("out"))
                        || (cc_toks.get(fwd + 1).copied() == Some("r")
                            && cc_toks.get(fwd + 2).copied() == Some("o")))
                    && cc_toks
                        .get(fwd + 3)
                        .is_some_and(|t| !CLAUSE_BREAKERS.contains(t));
                let is_pre_neg_intro = (matches!(tok, "negative" | "neg") && next_is_prep)
                    || (tok == "no"
                        && cc_toks.get(fwd + 1).copied() == Some("evidence")
                        && cc_toks.get(fwd + 2).copied() == Some("of"))
                    || (tok == "free" && cc_toks.get(fwd + 1).copied() == Some("of"))
                    || is_to_rule_out;
                let is_two_tok_breaker =
                    tok == "c" && cc_toks.get(fwd + 1).copied() == Some("o");
                let is_positive_marker = POSITIVE_POLARITY_MARKERS.contains(&tok)
                    && !(tok == "present" && fwd > post_start && cc_toks[fwd - 1] == "not");
                if CLAUSE_BREAKERS.contains(&tok)
                    || is_two_tok_breaker
                    || is_positive_marker
                    || matches!(tok, "denies" | "denied" | "deny" | "without")
                    || is_pre_neg_intro
                {
                    post_end = fwd;
                    break;
                }
            }
            let post_tokens = &cc_toks[post_start..post_end];
            let is_post_neg = POST_NEGATION_TRIGGERS.iter().any(|trig| {
                let trig_toks: Vec<&str> = trig.split_whitespace().collect();
                if trig_toks.is_empty() || post_tokens.len() < trig_toks.len() {
                    return false;
                }
                (0..=post_tokens.len() - trig_toks.len()).any(|i| {
                    post_tokens[i..i + trig_toks.len()] == trig_toks[..]
                        && !matches!(
                            cc_toks.get(post_start + i + trig_toks.len()).copied(),
                            Some("for" | "of")
                        )
                        && !(i > 0
                            && post_tokens[i - 1] == "not"
                            && matches!(
                                *trig,
                                "ruled out"
                                    | "rule out"
                                    | "excluded"
                                    | "resolved"
                                    | "absent"
                                    | "negative"
                                    | "neg"
                                    | "unlikely"
                            ))
                })
            });

            if is_pre_neg || is_post_neg {
                found_negated = true;
            } else {
                found_affirmative = true;
            }
        }
    }

    (found_affirmative, found_negated)
}

fn normalize_icd_list(icd_raw: &str) -> Vec<String> {
    icd_raw
        .split(|c: char| c == ';' || c == ',' || c.is_whitespace() || c == '|')
        .map(|s| s.trim().replace('.', "").to_ascii_uppercase())
        .filter(|s| !s.is_empty())
        .collect()
}

fn has_icd_prefix(codes: &[String], prefixes: &[&str]) -> Option<String> {
    for c in codes {
        for &p in prefixes {
            if c.starts_with(p) {
                return Some(format!("ICD:{c}"));
            }
        }
    }
    None
}

/// Recursive-descent boolean AST for NSSP ESSENCE CCDD queries supporting
/// `AND`, `OR`, `ANDNOT`, unary `NOT`, nested parentheses `(...)`, `,` (`OR` shorthand),
/// `^...^` substring/ICD-prefix terms, and word/token-bounded bare or quoted terms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EssenceQueryAst {
    Term(String),
    SubstringTerm(String),
    Not(Box<EssenceQueryAst>),
    And(Box<EssenceQueryAst>, Box<EssenceQueryAst>),
    Or(Box<EssenceQueryAst>, Box<EssenceQueryAst>),
    AndNot(Box<EssenceQueryAst>, Box<EssenceQueryAst>),
}

fn matches_word_bounded(haystack: &str, needle: &str) -> bool {
    let needle_toks: Vec<String> = needle
        .split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_ascii_lowercase())
        .collect();
    if needle_toks.is_empty() {
        return !needle.is_empty() && haystack.contains(needle);
    }
    let hay_toks: Vec<String> = haystack
        .split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_ascii_lowercase())
        .collect();
    if hay_toks.len() < needle_toks.len() {
        return false;
    }
    hay_toks
        .windows(needle_toks.len())
        .any(|w| w.iter().zip(needle_toks.iter()).all(|(a, b)| a == b))
}

impl EssenceQueryAst {
    pub fn evaluate(&self, padded_ccdd: &str) -> bool {
        match self {
            Self::Term(t) => !t.is_empty() && matches_word_bounded(padded_ccdd, t.as_str()),
            Self::SubstringTerm(t) => !t.is_empty() && padded_ccdd.contains(t.as_str()),
            Self::Not(inner) => !inner.evaluate(padded_ccdd),
            Self::And(lhs, rhs) => lhs.evaluate(padded_ccdd) && rhs.evaluate(padded_ccdd),
            Self::Or(lhs, rhs) => lhs.evaluate(padded_ccdd) || rhs.evaluate(padded_ccdd),
            Self::AndNot(lhs, rhs) => lhs.evaluate(padded_ccdd) && !rhs.evaluate(padded_ccdd),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum EssenceToken {
    LParen,
    RParen,
    And,
    Or,
    AndNot,
    Not,
    Term(String),
    SubstringTerm(String),
}

fn tokenize_essence_query(query_expr: &str) -> Result<Vec<EssenceToken>, EpiError> {
    let chars: Vec<char> = query_expr.chars().collect();
    let n = chars.len();
    let mut i = 0usize;
    let mut tokens: Vec<EssenceToken> = Vec::new();
    let mut bare_words: Vec<String> = Vec::new();

    let flush_bare = |words: &mut Vec<String>, out: &mut Vec<EssenceToken>| {
        if !words.is_empty() {
            out.push(EssenceToken::Term(words.join(" ")));
            words.clear();
        }
    };

    while i < n {
        let ch = chars[i];
        if ch.is_whitespace() {
            i += 1;
            continue;
        }
        if ch == '(' {
            flush_bare(&mut bare_words, &mut tokens);
            tokens.push(EssenceToken::LParen);
            i += 1;
            continue;
        }
        if ch == ')' {
            flush_bare(&mut bare_words, &mut tokens);
            tokens.push(EssenceToken::RParen);
            i += 1;
            continue;
        }
        if ch == ',' {
            flush_bare(&mut bare_words, &mut tokens);
            tokens.push(EssenceToken::Or);
            i += 1;
            continue;
        }
        if ch == '^' || ch == '"' || ch == '\'' {
            flush_bare(&mut bare_words, &mut tokens);
            let delim = ch;
            i += 1;
            let start = i;
            while i < n && chars[i] != delim {
                i += 1;
            }
            if i >= n {
                return Err(EpiError::Schema(format!(
                    "Unterminated {delim}...{delim} term in ESSENCE query"
                )));
            }
            let raw_inner: String = chars[start..i].iter().collect();
            i += 1;
            let cleaned = if delim == '^' {
                raw_inner.trim().replace('.', "").to_ascii_lowercase()
            } else {
                raw_inner
                    .split_whitespace()
                    .collect::<Vec<_>>()
                    .join(" ")
                    .to_ascii_lowercase()
            };
            if cleaned.is_empty() {
                return Err(EpiError::Schema(
                    "Empty delimited term in ESSENCE query".to_string(),
                ));
            }
            if delim == '^' {
                tokens.push(EssenceToken::SubstringTerm(cleaned));
            } else {
                tokens.push(EssenceToken::Term(cleaned));
            }
            continue;
        }

        let start = i;
        while i < n
            && !chars[i].is_whitespace()
            && !matches!(chars[i], '(' | ')' | ',' | '^' | '"' | '\'')
        {
            i += 1;
        }
        let word: String = chars[start..i].iter().collect();
        let lower = word.to_ascii_lowercase();
        match lower.as_str() {
            "and" => {
                flush_bare(&mut bare_words, &mut tokens);
                tokens.push(EssenceToken::And);
            }
            "or" => {
                flush_bare(&mut bare_words, &mut tokens);
                tokens.push(EssenceToken::Or);
            }
            "andnot" => {
                flush_bare(&mut bare_words, &mut tokens);
                tokens.push(EssenceToken::AndNot);
            }
            "not" => {
                flush_bare(&mut bare_words, &mut tokens);
                tokens.push(EssenceToken::Not);
            }
            _ => {
                bare_words.push(lower);
            }
        }
    }
    flush_bare(&mut bare_words, &mut tokens);
    Ok(tokens)
}

struct EssenceParser<'a> {
    tokens: &'a [EssenceToken],
    pos: usize,
}

impl<'a> EssenceParser<'a> {
    fn peek(&self) -> Option<&EssenceToken> {
        self.tokens.get(self.pos)
    }

    fn next(&mut self) -> Option<&EssenceToken> {
        let t = self.tokens.get(self.pos);
        if t.is_some() {
            self.pos += 1;
        }
        t
    }

    fn parse_or(&mut self) -> Result<EssenceQueryAst, EpiError> {
        let mut node = self.parse_and()?;
        while matches!(self.peek(), Some(EssenceToken::Or)) {
            self.next();
            let rhs = self.parse_and()?;
            node = EssenceQueryAst::Or(Box::new(node), Box::new(rhs));
        }
        Ok(node)
    }

    fn parse_and(&mut self) -> Result<EssenceQueryAst, EpiError> {
        let mut node = self.parse_not()?;
        loop {
            match self.peek() {
                Some(EssenceToken::And) => {
                    self.next();
                    let rhs = self.parse_not()?;
                    node = EssenceQueryAst::And(Box::new(node), Box::new(rhs));
                }
                Some(EssenceToken::AndNot) => {
                    self.next();
                    let rhs = self.parse_not()?;
                    node = EssenceQueryAst::AndNot(Box::new(node), Box::new(rhs));
                }
                _ => break,
            }
        }
        Ok(node)
    }

    fn parse_not(&mut self) -> Result<EssenceQueryAst, EpiError> {
        if matches!(self.peek(), Some(EssenceToken::Not)) {
            self.next();
            let inner = self.parse_not()?;
            return Ok(EssenceQueryAst::Not(Box::new(inner)));
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<EssenceQueryAst, EpiError> {
        match self.next().cloned() {
            Some(EssenceToken::LParen) => {
                let inner = self.parse_or()?;
                match self.next() {
                    Some(EssenceToken::RParen) => Ok(inner),
                    _ => Err(EpiError::Schema(
                        "Unclosed '(' in ESSENCE boolean query".to_string(),
                    )),
                }
            }
            Some(EssenceToken::Term(t)) => Ok(EssenceQueryAst::Term(t)),
            Some(EssenceToken::SubstringTerm(t)) => Ok(EssenceQueryAst::SubstringTerm(t)),
            Some(other) => Err(EpiError::Schema(format!(
                "Unexpected token {other:?} in ESSENCE boolean query"
            ))),
            None => Err(EpiError::Schema(
                "Unexpected end of ESSENCE boolean query".to_string(),
            )),
        }
    }
}

pub fn parse_essence_query(query_expr: &str) -> Result<EssenceQueryAst, EpiError> {
    let tokens = tokenize_essence_query(query_expr)?;
    if tokens.is_empty() {
        return Err(EpiError::Schema("Empty ESSENCE query".to_string()));
    }
    let mut parser = EssenceParser {
        tokens: &tokens,
        pos: 0,
    };
    let ast = parser.parse_or()?;
    if parser.pos != tokens.len() {
        return Err(EpiError::Schema(format!(
            "Unexpected trailing token {:?} in ESSENCE query",
            tokens[parser.pos]
        )));
    }
    Ok(ast)
}

/// Evaluate a custom NSSP ESSENCE boolean query supporting `AND`, `OR`, `ANDNOT`, `NOT`,
/// nested `(...)`, `,`, and `^term^` against concatenated normalized CC + ICD string.
pub fn evaluate_custom_essence_query(normalized_ccdd: &str, query_expr: &str) -> bool {
    let Ok(ast) = parse_essence_query(query_expr) else {
        return false;
    };
    let padded = format!(" {} ", normalized_ccdd.to_ascii_lowercase());
    ast.evaluate(&padded)
}

/// Classify a batch of ED triage encounters across all 6 CDC NSSP ESSENCE syndromes + optional custom query.
pub fn classify_ed_encounters_essence(
    encounters: &[EdTriageEncounter],
    custom_query_expr: &str,
    neg_window_tokens: usize,
) -> EssenceBatchClassification {
    let mut classified = Vec::with_capacity(encounters.len());
    let mut syn_counts: BTreeMap<EssenceSyndromeId, usize> = BTreeMap::new();
    let mut syn_neg_blocked: BTreeMap<EssenceSyndromeId, usize> = BTreeMap::new();
    for s in EssenceSyndromeId::ALL {
        syn_counts.insert(s, 0);
        syn_neg_blocked.insert(s, 0);
    }

    let mut any_matches = 0_usize;
    let mut total_neg_blocked = 0_usize;
    let mut custom_matches = 0_usize;

    for enc in encounters {
        let norm_cc = normalize_chief_complaint(&enc.chief_complaint_text);
        let icd_codes = normalize_icd_list(&enc.icd10_cm_codes);
        let ccdd_concat = format!("{} {}", norm_cc, icd_codes.join(" ").to_ascii_lowercase());

        let mut matched_syndromes = Vec::new();
        let mut matched_terms = Vec::new();
        let mut suppressed_terms = Vec::new();

        // Helper closure checking affirmative symptom list + global encounter exclusion phrases
        let mut check_terms = |syn: EssenceSyndromeId,
                               phrases: &[&str],
                               exclusions: &[&str]|
         -> Option<String> {
            if exclusions.iter().any(|ex| norm_cc.contains(ex)) {
                suppressed_terms.push(format!("{}:excluded_context", syn.code()));
                *syn_neg_blocked.entry(syn).or_insert(0) += 1;
                total_neg_blocked += 1;
                return None;
            }
            let mut hit_term: Option<String> = None;
            let mut only_neg = false;
            for &ph in phrases {
                let (aff, neg) = match_affirmative_or_negated(&norm_cc, ph, neg_window_tokens);
                if aff && hit_term.is_none() {
                    hit_term = Some(ph.to_string());
                } else if neg {
                    only_neg = true;
                    suppressed_terms.push(format!("negated({ph})"));
                }
            }
            if hit_term.is_none() && only_neg {
                *syn_neg_blocked.entry(syn).or_insert(0) += 1;
                total_neg_blocked += 1;
            }
            hit_term
        };

        // 1. CDC ILI: (fever AND (cough OR sore throat OR pharyngitis)) OR explicit influenza/ILI OR ICD J09-J11
        let ili_exclusions = &["immunization", "vaccine only", "flu shot"];
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["J09", "J10", "J11"]) {
            matched_syndromes.push(EssenceSyndromeId::IliRespiratory);
            matched_terms.push(icd_m);
        } else if let Some(direct_ili) = check_terms(
            EssenceSyndromeId::IliRespiratory,
            &["influenza", "flu like"],
            ili_exclusions,
        ) {
            matched_syndromes.push(EssenceSyndromeId::IliRespiratory);
            matched_terms.push(direct_ili);
        } else {
            let fev = check_terms(
                EssenceSyndromeId::IliRespiratory,
                &["fever", "febrile", "chills"],
                ili_exclusions,
            );
            let resp = check_terms(
                EssenceSyndromeId::IliRespiratory,
                &["cough", "sore throat", "pharyngitis"],
                ili_exclusions,
            );
            if let (Some(f_t), Some(r_t)) = (fev, resp) {
                matched_syndromes.push(EssenceSyndromeId::IliRespiratory);
                matched_terms.push(format!("{f_t}+{r_t}"));
            }
        }

        // 2. COVID-Like Illness (CLI): U071/J1282 OR (cough/shortness of breath/anosmia + fever/covid)
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["U071", "J1282"]) {
            matched_syndromes.push(EssenceSyndromeId::CovidLikeIllness);
            matched_terms.push(icd_m);
        } else if let Some(cli_hit) = check_terms(
            EssenceSyndromeId::CovidLikeIllness,
            &["covid", "sars cov", "loss of taste", "loss of smell", "anosmia"],
            &["covid test for travel", "vaccine"],
        ) {
            matched_syndromes.push(EssenceSyndromeId::CovidLikeIllness);
            matched_terms.push(cli_hit);
        } else {
            let dys = match_affirmative_or_negated(&norm_cc, "shortness of breath", neg_window_tokens).0
                || match_affirmative_or_negated(&norm_cc, "hypoxia", neg_window_tokens).0;
            let fev = match_affirmative_or_negated(&norm_cc, "fever", neg_window_tokens).0;
            if dys && fev {
                matched_syndromes.push(EssenceSyndromeId::CovidLikeIllness);
                matched_terms.push("fever+shortness_of_breath".to_string());
            }
        }

        // 3. Acute GI / Enteric: A08/A09/R112 OR vomiting/diarrhea/gastroenteritis (excluding pregnancy/chemo)
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["A08", "A09", "R197"]) {
            matched_syndromes.push(EssenceSyndromeId::AcuteGiEnteric);
            matched_terms.push(icd_m);
        } else if let Some(gi_hit) = check_terms(
            EssenceSyndromeId::AcuteGiEnteric,
            &["diarrhea", "vomiting", "gastroenteritis", "norovirus", "watery stool"],
            &["chemotherapy", "morning sickness", "hyperemesis gravidarum"],
        ) {
            matched_syndromes.push(EssenceSyndromeId::AcuteGiEnteric);
            matched_terms.push(gi_hit);
        }

        // 4. Suspected Unintentional Overdose (CDC DOSE): T40/T42/T43 (non-intentional) OR narcan/naloxone/overdose
        let od_exclusions = &["self harm", "suicide attempt", "intentional overdose", "detox clearance"];
        let mut od_intentional_excluded = false;
        if !od_exclusions.iter().any(|ex| norm_cc.contains(ex)) {
            if let Some(icd_m) = has_icd_prefix(&icd_codes, &["T400", "T401", "T402", "T404", "T405", "T436"]) {
                matched_syndromes.push(EssenceSyndromeId::SuspectedOverdose);
                matched_terms.push(icd_m);
            } else if let Some(od_hit) = check_terms(
                EssenceSyndromeId::SuspectedOverdose,
                &["overdose", "naloxone", "narcan", "fentanyl", "heroin", "pinpoint pupils"],
                od_exclusions,
            ) {
                matched_syndromes.push(EssenceSyndromeId::SuspectedOverdose);
                matched_terms.push(od_hit);
            }
        } else {
            od_intentional_excluded = true;
        }

        // 5. Heat-Related Illness (HRI): T67/X30 OR heat exhaustion / heat stroke / sunstroke (excluding "heat pack"!)
        let hri_exclusions = &["heat pack", "heating pad", "chili heat", "prickly heat"];
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["T67", "X30"]) {
            matched_syndromes.push(EssenceSyndromeId::HeatRelatedIllness);
            matched_terms.push(icd_m);
        } else if let Some(hri_hit) = check_terms(
            EssenceSyndromeId::HeatRelatedIllness,
            &[
                "heat exhaustion",
                "heat stroke",
                "sunstroke",
                "heat cramps",
                "working outside in heat",
                "overheated outdoors",
            ],
            hri_exclusions,
        ) {
            matched_syndromes.push(EssenceSyndromeId::HeatRelatedIllness);
            matched_terms.push(hri_hit);
        }

        // 6. Fever + Rash / Exanthem: B05/B06 OR (fever + rash / maculopapular / koplik)
        if let Some(icd_m) = has_icd_prefix(&icd_codes, &["B05", "B06"]) {
            matched_syndromes.push(EssenceSyndromeId::FeverRashExanthem);
            matched_terms.push(icd_m);
        } else {
            let rash = check_terms(
                EssenceSyndromeId::FeverRashExanthem,
                &["rash", "maculopapular", "exanthem", "koplik"],
                &["diaper rash", "poison oak", "ivy rash"],
            );
            let fev = match_affirmative_or_negated(&norm_cc, "fever", neg_window_tokens).0;
            if let Some(r_t) = rash {
                if fev || r_t == "koplik" {
                    matched_syndromes.push(EssenceSyndromeId::FeverRashExanthem);
                    matched_terms.push(format!("fever+{r_t}"));
                }
            }
        }

        if od_intentional_excluded {
            suppressed_terms.push("DOSE_OVERDOSE_V2:intentional_excluded".to_string());
            *syn_neg_blocked.entry(EssenceSyndromeId::SuspectedOverdose).or_insert(0) += 1;
            total_neg_blocked += 1;
        }

        for &syn in &matched_syndromes {
            *syn_counts.entry(syn).or_insert(0) += 1;
        }
        if !matched_syndromes.is_empty() {
            any_matches += 1;
        }

        let custom_hit = evaluate_custom_essence_query(&ccdd_concat, custom_query_expr);
        if custom_hit {
            custom_matches += 1;
        }

        classified.push(ClassifiedEdEncounter {
            visit_id: enc.visit_id.clone(),
            visit_iso_date: enc.visit_iso_date.clone(),
            facility_id: enc.facility_id.clone(),
            normalized_cc_text: norm_cc,
            matched_syndromes,
            matched_inclusion_terms: matched_terms,
            suppressed_negated_terms: suppressed_terms,
            custom_query_hit: custom_hit,
        });
    }

    let denom = (encounters.len().max(1)) as f64;
    let syndrome_summaries = EssenceSyndromeId::ALL
        .iter()
        .map(|&s| {
            let c = syn_counts.get(&s).copied().unwrap_or(0);
            SyndromeSummaryRow {
                syndrome: s,
                code: s.code().to_string(),
                label: s.label().to_string(),
                matched_visits: c,
                pct_of_ed_visits: (c as f64) * 100.0 / denom,
                negated_false_positives_blocked: syn_neg_blocked.get(&s).copied().unwrap_or(0),
            }
        })
        .collect();

    EssenceBatchClassification {
        total_encounters: encounters.len(),
        encounters_matching_any_syndrome: any_matches,
        total_negations_blocked: total_neg_blocked,
        custom_query_expr: custom_query_expr.to_string(),
        custom_query_matches: custom_matches,
        syndrome_summaries,
        classified_encounters: classified,
    }
}

/// Aggregate classified ED visits for a chosen syndrome into an MMWR-week `TimeSeriesDataset` for SQLite & EARS/Farrington.
pub fn build_syndrome_weekly_dataset(
    batch: &EssenceBatchClassification,
    syndrome: EssenceSyndromeId,
) -> TimeSeriesDataset {
    let mut week_buckets: BTreeMap<(i32, u8), (Option<epi_core::MmwrWeek>, usize)> = BTreeMap::new();
    for enc in &batch.classified_encounters {
        let Some(mmwr) = parse_flexible_date_to_mmwr(&enc.visit_iso_date) else {
            continue;
        };
        let entry = week_buckets
            .entry((mmwr.year, mmwr.week))
            .or_insert((Some(mmwr), 0));
        if enc.matched_syndromes.contains(&syndrome) {
            entry.1 += 1;
        }
    }

    let points: Vec<TimeSeriesPoint> = week_buckets
        .into_iter()
        .map(|((yr, wk), (mmwr_opt, cnt))| TimeSeriesPoint {
            date_label: format!("{yr}-W{wk:02}"),
            epiweek: mmwr_opt,
            count: cnt as f64,
            suppressed: false,
        })
        .collect();

    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for p in &points {
        for b in p.date_label.bytes() {
            hash ^= u64::from(b);
            hash = hash.wrapping_mul(0x0100_0000_01b3);
        }
        hash ^= p.count.to_bits();
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }

    TimeSeriesDataset {
        vintage_id: format!("nssp_essence_{}_v1", syndrome.code().to_ascii_lowercase()),
        title: format!("NSSP ESSENCE Weekly ED Visits: {}", syndrome.label()),
        source_kind: SourceKind::LocalCsvTsv,
        source_url: "local://nssp-essence-ccdd-classifier".to_string(),
        ingested_at_utc: "2026-11-13T00:00:00Z".to_string(),
        default_si_preset: "influenza".to_string(),
        checksum_hex: format!("{:08x}", (hash & 0xffff_ffff) as u32),
        points,
    }
}

/// Parse ED encounters CSV (`visit_id,visit_iso_date,facility_id,patient_age,chief_complaint_text,icd10_cm_codes`).
pub fn parse_ed_encounters_csv(csv_text: &str) -> Result<Vec<EdTriageEncounter>, EpiError> {
    let mut rows = Vec::new();
    for (idx, raw_line) in csv_text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        if idx == 0 && cols[0].eq_ignore_ascii_case("visit_id") {
            continue;
        }
        if cols.len() < 5 {
            return Err(EpiError::Schema(format!(
                "ED encounters CSV line {} requires >=5 columns",
                idx + 1
            )));
        }
        rows.push(EdTriageEncounter {
            visit_id: cols[0].to_string(),
            visit_iso_date: cols[1].to_string(),
            facility_id: cols[2].to_string(),
            patient_age: cols[3].parse().unwrap_or(35),
            chief_complaint_text: cols[4].to_string(),
            icd10_cm_codes: cols.get(5).copied().unwrap_or("").to_string(),
        });
    }
    if rows.is_empty() {
        return Err(EpiError::Schema(
            "No valid ED encounters found in CSV".to_string(),
        ));
    }
    Ok(rows)
}

/// Built-in 8-week 68-encounter ED syndromic triage dataset with clinical negations, exclusion idioms, and an ILI/GI wave.
pub fn builtin_ed_chief_complaint_encounters() -> Vec<EdTriageEncounter> {
    let mut out = Vec::new();
    let week_dates = [
        "2026-09-21",
        "2026-09-28",
        "2026-10-05",
        "2026-10-12",
        "2026-10-19",
        "2026-10-26",
        "2026-11-02",
        "2026-11-09",
    ];

    // Tricky negation & exclusion benchmark encounters (must NOT trigger false-positive syndrome counts)
    let tricky_negations: &[(&str, &str, &str)] = &[
        (
            "2026-09-22",
            "Lac right thumb; pt denies fever cough or sob",
            "S61.011A",
        ),
        (
            "2026-09-29",
            "Low back strain after lifting; applied heat pack at home without relief",
            "M54.50",
        ),
        (
            "2026-10-06",
            "Ankle inversion injury; no fever no chills negative for vomiting or diarrhea",
            "S93.401A",
        ),
        (
            "2026-10-13",
            "Here for routine flu shot immunization; denies cough or sore throat",
            "Z23",
        ),
        (
            "2026-10-20",
            "Intentional overdose self harm evaluation by psychiatry",
            "T43.202A",
        ),
        (
            "2026-10-27",
            "8 wks pregnant morning sickness hyperemesis gravidarum; denies fever",
            "O21.0",
        ),
        (
            "2026-11-03",
            "Contact poison oak rash on forearms after trail work; no fever",
            "L23.7",
        ),
        (
            "2026-11-10",
            "Rule out DVT left calf; denies sob denies chest pain no fever",
            "I82.402",
        ),
    ];

    for (idx, &(dt, cc, icd)) in tricky_negations.iter().enumerate() {
        out.push(EdTriageEncounter {
            visit_id: format!("ED-NEG-{:03}", idx + 1),
            visit_iso_date: dt.to_string(),
            facility_id: "OR-GOODSAM-ED".to_string(),
            patient_age: 29 + (idx as u8) * 4,
            chief_complaint_text: cc.to_string(),
            icd10_cm_codes: icd.to_string(),
        });
    }

    // Weekly background + escalating Weeks 6-8 ILI and Acute GI outbreak signal
    let ili_per_week = [1_usize, 1, 2, 2, 3, 5, 7, 10];
    let gi_per_week = [1_usize, 1, 1, 2, 2, 4, 5, 6];
    let mut seq = 1_usize;

    for (w_idx, &dt) in week_dates.iter().enumerate() {
        for k in 0..ili_per_week[w_idx] {
            let (cc, icd) = if k % 3 == 0 {
                ("High fever chills myalgias and productive cough x2 days", "J10.1")
            } else if k % 3 == 1 {
                ("Tmax 102.4 sore throat cough; denies vomiting", "R50.9;J02.9")
            } else {
                ("Fever sob and loss of smell after household covid exposure", "U07.1")
            };
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-GOODSAM-ED".to_string(),
                patient_age: 19 + ((seq * 7) % 62) as u8,
                chief_complaint_text: cc.to_string(),
                icd10_cm_codes: icd.to_string(),
            });
            seq += 1;
        }

        for k in 0..gi_per_week[w_idx] {
            let (cc, icd) = if k % 2 == 0 {
                ("Acute n/v/d watery stool after banquet; denies cough", "A08.4")
            } else {
                ("Intractable emesis and diarrhea x18h", "A09")
            };
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-RIVERBEND-ED".to_string(),
                patient_age: 12 + ((seq * 5) % 58) as u8,
                chief_complaint_text: cc.to_string(),
                icd10_cm_codes: icd.to_string(),
            });
            seq += 1;
        }

        // Add sentinel overdose, heat, and fever-rash alerts in selected weeks
        if w_idx % 2 == 1 {
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-RIVERBEND-ED".to_string(),
                patient_age: 34,
                chief_complaint_text:
                    "Found unresponsive with pinpoint pupils; EMS gave narcan 4mg intranasal"
                        .to_string(),
                icd10_cm_codes: "T40.411A".to_string(),
            });
            seq += 1;
        }
        if w_idx == 0 {
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-GOODSAM-ED".to_string(),
                patient_age: 46,
                chief_complaint_text:
                    "Dizziness and heat exhaustion after working outside in heat harvesting"
                        .to_string(),
                icd10_cm_codes: "T67.5XXA".to_string(),
            });
            seq += 1;
        }
        if w_idx >= 6 {
            out.push(EdTriageEncounter {
                visit_id: format!("ED-VISIT-{seq:03}"),
                visit_iso_date: dt.to_string(),
                facility_id: "OR-GOODSAM-ED".to_string(),
                patient_age: 7,
                chief_complaint_text:
                    "Unvaccinated child with high fever coryza and maculopapular rash spreading from face"
                        .to_string(),
                icd10_cm_codes: "B05.9".to_string(),
            });
            seq += 1;
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_negex_blocks_denied_symptoms_and_heat_pack_idiom() {
        let encounters = builtin_ed_chief_complaint_encounters();
        let batch = classify_ed_encounters_essence(
            &encounters,
            "(^fever^ OR ^cough^) ANDNOT (^immunization^)",
            5,
        );

        // First 8 encounters are all negated or excluded idioms -> zero syndrome matches!
        for neg_enc in batch.classified_encounters.iter().take(8) {
            assert!(
                neg_enc.matched_syndromes.is_empty(),
                "Negated/excluded encounter {} falsely matched {:?}",
                neg_enc.visit_id,
                neg_enc.matched_syndromes
            );
        }
        assert!(batch.total_negations_blocked >= 6);

        let ds = build_syndrome_weekly_dataset(&batch, EssenceSyndromeId::IliRespiratory);
        assert_eq!(ds.points.len(), 8);
        assert!(ds.points.last().unwrap().count > ds.points.first().unwrap().count);
    }

    #[test]
    fn test_post_negations_and_clause_isolation() {
        let norm1 = normalize_chief_complaint("Fever ruled out; cough present");
        let (fev_aff, fev_neg) = match_affirmative_or_negated(&norm1, "fever", 5);
        assert!(!fev_aff, "fever ruled out must not be affirmative");
        assert!(fev_neg, "fever ruled out must record negation");
        let (cough_aff, cough_neg) = match_affirmative_or_negated(&norm1, "cough", 5);
        assert!(cough_aff, "cough in next clause must stay affirmative");
        assert!(!cough_neg);

        let norm2 = normalize_chief_complaint("covid:negative; influenza not detected; fever positive");
        let (cov_aff, cov_neg) = match_affirmative_or_negated(&norm2, "covid", 5);
        assert!(!cov_aff && cov_neg, "colon-delimited post-negation must block covid");
        let (flu_aff, flu_neg) = match_affirmative_or_negated(&norm2, "influenza", 5);
        assert!(!flu_aff && flu_neg, "not detected post-negation must block influenza");
        let (fev2_aff, _) = match_affirmative_or_negated(&norm2, "fever", 5);
        assert!(fev2_aff, "fever positive after semicolon must remain affirmative");

        // Post-negation without semicolon must not bleed into following positive symptom
        let norm3 = normalize_chief_complaint("cough not detected fever positive");
        let (c3_aff, c3_neg) = match_affirmative_or_negated(&norm3, "cough", 5);
        let (f3_aff, f3_neg) = match_affirmative_or_negated(&norm3, "fever", 5);
        assert!(!c3_aff && c3_neg);
        assert!(f3_aff && !f3_neg);
    }

    #[test]
    fn test_essence_boolean_ast_precedence_nesting_and_word_boundaries() {
        // Word-boundary safety: "or" inside "sore throat" and "norovirus" must NOT split
        let cc1 = "sore throat and fever j101";
        assert!(evaluate_custom_essence_query(
            cc1,
            "(sore throat OR norovirus) AND ^J10^"
        ));
        assert!(!evaluate_custom_essence_query(
            "sore throat only r070",
            "(sore throat OR norovirus) AND ^J10^"
        ));
        assert!(evaluate_custom_essence_query(
            "acute norovirus gastroenteritis",
            "sore throat OR norovirus"
        ));

        // Nested precedence: AND binds tighter than OR; NOT and ANDNOT exclude correctly
        let q = "((^fever^ AND (^cough^ OR ^dyspnea^)) OR ^anosmia^) ANDNOT (^vaccine^ OR NOT ^cough^)";
        assert!(evaluate_custom_essence_query("fever and cough u071", q));
        assert!(!evaluate_custom_essence_query("fever and cough vaccine", q));
        assert!(!evaluate_custom_essence_query("anosmia without dyspnea", q));
        assert!(evaluate_custom_essence_query("anosmia", "^anosmia^ OR ^fever^ AND ^cough^"));
        assert!(!evaluate_custom_essence_query("fever", "^anosmia^ OR ^fever^ AND ^cough^"));

        // Malformed queries return false safely without panicking
        assert!(!evaluate_custom_essence_query(cc1, ""));
        assert!(!evaluate_custom_essence_query(cc1, "(^fever^ AND"));
        assert!(!evaluate_custom_essence_query(cc1, "^unclosed"));

        // Bare/quoted terms must respect word/token boundaries while ^...^ keeps substring matching
        assert!(!evaluate_custom_essence_query("acid reflux", "flu"));
        assert!(!evaluate_custom_essence_query("motor vehicle crash", "rash"));
        assert!(!evaluate_custom_essence_query("blood in stool", "od"));
        assert!(!evaluate_custom_essence_query("patient has pain", "ha"));
        assert!(!evaluate_custom_essence_query("acid reflux", "\"flu\""));
        assert!(evaluate_custom_essence_query("acute flu symptoms", "flu"));
        assert!(evaluate_custom_essence_query("maculopapular rash", "rash"));
        assert!(evaluate_custom_essence_query("high fever j101", "^j10^"));
    }

    #[test]
    fn test_negex_challenger_regressions_and_normalization() {
        // 1. "no evidence of <symptom>" must negate the target symptom
        let norm_measles = normalize_chief_complaint("no evidence of measles");
        assert_eq!(
            match_affirmative_or_negated(&norm_measles, "measles", 5),
            (false, true),
            "no evidence of measles must be negated"
        );
        let norm_measles_fev = normalize_chief_complaint("no evidence of measles; fever present");
        assert_eq!(
            match_affirmative_or_negated(&norm_measles_fev, "measles", 5),
            (false, true)
        );
        assert_eq!(
            match_affirmative_or_negated(&norm_measles_fev, "fever", 5),
            (true, false)
        );
        // Post-negation form "measles no evidence" still negates
        let norm_post_ev = normalize_chief_complaint("measles no evidence; cough present");
        assert_eq!(
            match_affirmative_or_negated(&norm_post_ev, "measles", 5),
            (false, true)
        );
        assert_eq!(
            match_affirmative_or_negated(&norm_post_ev, "cough", 5),
            (true, false)
        );

        // 2. Punctuation-free post-negation bleed & "negative for" false post-negation
        let n_covid_mild = normalize_chief_complaint("covid negative mild fever");
        assert_eq!(
            match_affirmative_or_negated(&n_covid_mild, "covid", 5),
            (false, true)
        );
        assert_eq!(
            match_affirmative_or_negated(&n_covid_mild, "fever", 5),
            (true, false)
        );

        for text in [
            "fever absent cough present",
            "fever resolved cough present",
            "fever rule out cough present",
        ] {
            let n = normalize_chief_complaint(text);
            assert_eq!(
                match_affirmative_or_negated(&n, "fever", 5),
                (false, true),
                "expected fever negated in {text:?} (norm={n:?})"
            );
            assert_eq!(
                match_affirmative_or_negated(&n, "cough", 5),
                (true, false),
                "expected cough affirmative in {text:?} (norm={n:?})"
            );
        }

        let n_neg_for = normalize_chief_complaint("has cough, negative for fever");
        assert_eq!(
            match_affirmative_or_negated(&n_neg_for, "cough", 5),
            (true, false),
            "cough must remain affirmative before 'negative for fever'"
        );
        assert_eq!(
            match_affirmative_or_negated(&n_neg_for, "fever", 5),
            (false, true),
            "fever must be negated by 'negative for'"
        );

        // 3. Normalization: decimals, parens/brackets, contractions
        let n_dec = normalize_chief_complaint("Tmax 102.4 ruled out");
        assert_eq!(n_dec, "fever 102.4 ruled out");
        assert_eq!(
            match_affirmative_or_negated(&n_dec, "fever", 5),
            (false, true)
        );

        let n_paren = normalize_chief_complaint("fever(ruled out)");
        assert_eq!(n_paren, "fever ruled out");
        assert_eq!(
            match_affirmative_or_negated(&n_paren, "fever", 5),
            (false, true)
        );

        let n_contr = normalize_chief_complaint("don't have fever");
        assert_eq!(
            match_affirmative_or_negated(&n_contr, "fever", 5),
            (false, true)
        );
    }

    #[test]
    fn test_comma_pre_negation_boundaries_and_co_clause_breaker() {
        // 1. "has cough, rule out covid" -> cough affirmative, covid negated, no CovidLikeIllness
        let enc_ro = EdTriageEncounter {
            visit_id: "ADV-1".to_string(),
            visit_iso_date: "2026-11-10".to_string(),
            facility_id: "OR-GOODSAM-ED".to_string(),
            patient_age: 42,
            chief_complaint_text: "has cough, rule out covid".to_string(),
            icd10_cm_codes: String::new(),
        };
        let norm_ro = normalize_chief_complaint(&enc_ro.chief_complaint_text);
        assert_eq!(match_affirmative_or_negated(&norm_ro, "cough", 5), (true, false));
        assert_eq!(match_affirmative_or_negated(&norm_ro, "covid", 5), (false, true));
        let batch_ro = classify_ed_encounters_essence(&[enc_ro], "cough", 5);
        assert!(
            !batch_ro.classified_encounters[0]
                .matched_syndromes
                .contains(&EssenceSyndromeId::CovidLikeIllness),
            "ruled-out covid must not trigger CovidLikeIllness"
        );

        // 2. "fever, r/o covid" -> fever affirmative, covid negated, no CovidLikeIllness
        let enc_sl = EdTriageEncounter {
            visit_id: "ADV-2".to_string(),
            visit_iso_date: "2026-11-10".to_string(),
            facility_id: "OR-GOODSAM-ED".to_string(),
            patient_age: 35,
            chief_complaint_text: "fever, r/o covid".to_string(),
            icd10_cm_codes: String::new(),
        };
        let norm_sl = normalize_chief_complaint(&enc_sl.chief_complaint_text);
        assert_eq!(match_affirmative_or_negated(&norm_sl, "fever", 5), (true, false));
        assert_eq!(match_affirmative_or_negated(&norm_sl, "covid", 5), (false, true));
        let batch_sl = classify_ed_encounters_essence(&[enc_sl], "fever", 5);
        assert!(
            !batch_sl.classified_encounters[0]
                .matched_syndromes
                .contains(&EssenceSyndromeId::CovidLikeIllness),
            "r/o covid must not trigger CovidLikeIllness"
        );

        // 3. "fever, absent cough" -> fever affirmative, cough negated
        let norm_abs = normalize_chief_complaint("fever, absent cough");
        assert_eq!(match_affirmative_or_negated(&norm_abs, "fever", 5), (true, false));
        assert_eq!(match_affirmative_or_negated(&norm_abs, "cough", 5), (false, true));

        // 4. "denies fever c/o cough" (and un-collapsed "denies fever c o cough") -> fever negated, cough affirmative
        let norm_co = normalize_chief_complaint("denies fever c/o cough");
        assert_eq!(match_affirmative_or_negated(&norm_co, "fever", 5), (false, true));
        assert_eq!(match_affirmative_or_negated(&norm_co, "cough", 5), (true, false));
        assert_eq!(
            match_affirmative_or_negated("denies fever c o cough", "fever", 5),
            (false, true)
        );
        assert_eq!(
            match_affirmative_or_negated("denies fever c o cough", "cough", 5),
            (true, false)
        );

        // 5. Coordinate negated list "no fever, chills, cough" -> all three negated
        let norm_list = normalize_chief_complaint("no fever, chills, cough");
        assert_eq!(match_affirmative_or_negated(&norm_list, "fever", 5), (false, true));
        assert_eq!(match_affirmative_or_negated(&norm_list, "chills", 5), (false, true));
        assert_eq!(match_affirmative_or_negated(&norm_list, "cough", 5), (false, true));

        // 6. Edge cases: "fever, not ruled out" stays affirmative; "fever, ruled out" and "covid, negative" stay negated
        let norm_nro = normalize_chief_complaint("fever, not ruled out");
        assert_eq!(match_affirmative_or_negated(&norm_nro, "fever", 5), (true, false));
        let norm_ro_post = normalize_chief_complaint("fever, ruled out");
        assert_eq!(match_affirmative_or_negated(&norm_ro_post, "fever", 5), (false, true));
        let norm_neg_post = normalize_chief_complaint("covid, negative");
        assert_eq!(match_affirmative_or_negated(&norm_neg_post, "covid", 5), (false, true));
    }

    #[test]
    fn test_not_present_post_negation_and_diag_edge_cases() {
        // Primary blocker: "fever, not present" and "fever not present" must negate fever
        let norm_comma_np = normalize_chief_complaint("fever, not present");
        assert_eq!(
            match_affirmative_or_negated(&norm_comma_np, "fever", 5),
            (false, true),
            "fever, not present must be negated"
        );
        let norm_plain_np = normalize_chief_complaint("fever not present");
        assert_eq!(
            match_affirmative_or_negated(&norm_plain_np, "fever", 5),
            (false, true),
            "fever not present must be negated"
        );
        let norm_np_then_pos = normalize_chief_complaint("fever not present, cough present");
        assert_eq!(match_affirmative_or_negated(&norm_np_then_pos, "fever", 5), (false, true));
        assert_eq!(match_affirmative_or_negated(&norm_np_then_pos, "cough", 5), (true, false));

        // Positive "present" and preceding post-negation + "present" must remain affirmative
        let norm_pos_pres = normalize_chief_complaint("fever present");
        assert_eq!(
            match_affirmative_or_negated(&norm_pos_pres, "fever", 5),
            (true, false),
            "fever present must remain affirmative"
        );
        let norm_nd_fp = normalize_chief_complaint("cough not detected fever present");
        assert_eq!(match_affirmative_or_negated(&norm_nd_fp, "cough", 5), (false, true));
        assert_eq!(match_affirmative_or_negated(&norm_nd_fp, "fever", 5), (true, false));

        // DIAG-1: "has cough, to rule out covid" (and unpunctuated / r/o variants)
        for text in [
            "has cough, to rule out covid",
            "has cough to rule out covid",
            "has cough, to r/o covid",
        ] {
            let n = normalize_chief_complaint(text);
            assert_eq!(
                match_affirmative_or_negated(&n, "cough", 5),
                (true, false),
                "cough must stay affirmative in {text:?} (norm={n:?})"
            );
            assert_eq!(
                match_affirmative_or_negated(&n, "covid", 5),
                (false, true),
                "covid must be negated in {text:?} (norm={n:?})"
            );
        }

        // DIAG-2: "fever, negative for influenza, positive for covid" (and negated positive marker)
        let norm_diag2 = normalize_chief_complaint("fever, negative for influenza, positive for covid");
        assert_eq!(match_affirmative_or_negated(&norm_diag2, "fever", 5), (true, false));
        assert_eq!(match_affirmative_or_negated(&norm_diag2, "influenza", 5), (false, true));
        assert_eq!(match_affirmative_or_negated(&norm_diag2, "covid", 5), (true, false));

        let norm_not_pos = normalize_chief_complaint("not positive for covid");
        assert_eq!(match_affirmative_or_negated(&norm_not_pos, "covid", 5), (false, true));
    }
}
