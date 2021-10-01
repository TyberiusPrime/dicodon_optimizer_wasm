#![allow(dead_code)]
mod utils;

//use seq_io::fasta::{Reader, Record};
use itertools::Itertools;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use wasm_bindgen::prelude::*;

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

#[wasm_bindgen]
extern "C" {
    fn alert(s: &str);
}

extern crate web_sys;

// A macro to provide `println!(..)`-style syntax for `console.log` logging.
macro_rules! log {
    ( $( $t:tt )* ) => {
        web_sys::console::log_1(&format!( $( $t )* ).into());
    }
}

#[wasm_bindgen]
pub fn greet(name: &str) {
    alert(&format!("Hello {}! How are you?", name));
}

#[wasm_bindgen]
pub fn start() {
    utils::set_panic_hook();
}

#[wasm_bindgen]
pub fn fasta_to_score_table(fasta: &str) -> String {
    if !fasta.starts_with(">") {
        return "File was not fasta. Did not start with >".to_string();
    }
    let counts = fasta_to_count_table(fasta);
    let normed = norm_count_table_by_di_aa(counts);
    //format!("{:?}",
    pretty_dicodon_table(normed)
}

/// count dicodons, in 6*3 bit encoding,
/// internally using a 1MB vector
/// convert to lower case string for output
fn fasta_to_count_table(fasta: &str) -> HashMap<String, u32> {
    let mut counter: Vec<u32> = vec![0; 1048576];
    'outer: for subsection in fasta.split("\n>") {
        let subbytes = subsection.as_bytes();
        let it = subbytes.iter();
        let mut it = it.skip_while(|&&x| x != b'\n');
        it.next(); //eat \n
        let mut dicodon: u32 = 0;
        for _ in 0..6 {
            let x = it.next();
            match x {
                Some(x) => {
                    dicodon = dicodon << 3 | encode_base(*x);
                }
                None => {
                    log!("too short");
                    continue 'outer; //too short...
                }
            }
        }
        for base in it {
            counter[dicodon as usize] += 1;
            dicodon = (dicodon << 3 | encode_base(*base)) & 0b111111111111111111;
            // 6 bases..
        }
    }
    let mut dicodon_counter = HashMap::<String, u32>::new();
    for (enc_dicodon, count) in counter.into_iter().enumerate() {
        if count > 0 {
            let sdicodon = decode_dicodon(enc_dicodon as u32);
            if !sdicodon.contains('N') {
                dicodon_counter.insert(sdicodon, count);
            }
        }
    }
    dicodon_counter
}

fn norm_count_table_by_di_aa(counts: HashMap<String, u32>) -> HashMap<String, (f64, u32)> {
    let mut res = HashMap::new();
    for (_diaa, dicodons) in dicodons_by_amino_acids().iter() {
        let total: u32 = dicodons.iter().map(|x| counts.get(x).unwrap_or(&0)).sum();
        for d in dicodons.iter() {
            let dicodon_count = counts.get(d).unwrap_or(&0);
            res.insert(
                d.to_string(),
                ((*dicodon_count as f64) / (total as f64), *dicodon_count),
            );
        }
    }

    res
}

fn pretty_dicodon_table(input: HashMap<String, (f64, u32)>) -> String {
    let mut res = Vec::new();
    for (diaa, dicodons) in dicodons_by_amino_acids().iter() {
        for dicodon in dicodons.iter() {
            let e = input.get(dicodon).unwrap_or(&(0., 0));
            res.push(format!("{}\t{}\t{}\t{}", diaa, dicodon, e.0, e.1));
        }
    }
    res.sort();
    format!("AA\tDicodon\tRel.Freq\tCount\n{}", res.join("\n"))
}

fn encode_base(b: u8) -> u32 {
    match b {
        b'A' => 0b001,
        b'a' => 0b001,
        b'C' => 0b010,
        b'c' => 0b010,
        b'G' => 0b011,
        b'g' => 0b011,
        b'T' => 0b100,
        b't' => 0b100,
        _ => 0b111,
    }
}

fn decode_base(eb: u32) -> u8 {
    match eb {
        0b001 => b'A',
        0b010 => b'C',
        0b011 => b'G',
        0b100 => b'T',
        _ => b'N',
    }
}

fn decode_dicodon(dicodon: u32) -> String {
    let f1 = decode_base(dicodon & 0b111);
    let f2 = decode_base(dicodon >> 3 & 0b111);
    let f3 = decode_base(dicodon >> 6 & 0b111);
    let f4 = decode_base(dicodon >> 9 & 0b111);
    let f5 = decode_base(dicodon >> 12 & 0b111);
    let f6 = decode_base(dicodon >> 15 & 0b111);
    let mut r: Vec<u8> = Vec::new();
    r.push(f6);
    r.push(f5);
    r.push(f4);
    r.push(f3);
    r.push(f2);
    r.push(f1);
    std::str::from_utf8(&r).unwrap().to_string()
}

macro_rules! collection {
    // map-like
    ($($k:expr => $v:expr),* $(,)?) => {{
        use std::iter::{Iterator, IntoIterator};
        Iterator::collect(IntoIterator::into_iter([$(($k.to_owned(), $v.to_owned()),)*]))
    }};
}

fn codon_to_aa_map() -> HashMap<String, String> {
    collection! {
    "AAA" => "K",
    "AAC" => "N",
    "AAG" => "K",
    "AAT" => "N",
    "ACA" => "T",
    "ACC" => "T",
    "ACG" => "T",
    "ACT" => "T",
    "AGA" => "R",
    "AGC" => "S",
    "AGG" => "R",
    "AGT" => "S",
    "ATA" => "I",
    "ATC" => "I",
    "ATG" => "M",
    "ATT" => "I",
    "ATT" => "I",
    "CAA" => "Q",
    "CAC" => "H",
    "CAG" => "Q",
    "CAT" => "H",
    "CCA" => "P",
    "CCC" => "P",
    "CCG" => "P",
    "CCT" => "P",
    "CGA" => "R",
    "CGC" => "R",
    "CGG" => "R",
    "CGT" => "R",
    "CTA" => "L",
    "CTC" => "L",
    "CTG" => "L",
    "CTT" => "L",
    "GAA" => "E",
    "GAC" => "D",
    "GAG" => "E",
    "GAT" => "D",
    "GCA" => "A",
    "GCC" => "A",
    "GCG" => "A",
    "GCT" => "A",
    "GGA" => "G",
    "GGC" => "G",
    "GGG" => "G",
    "GGT" => "G",
    "GTA" => "V",
    "GTC" => "V",
    "GTG" => "V",
    "GTT" => "V",
    "TAA" => "*",
    "TAC" => "Y",
    "TAG" => "*",
    "TAT" => "Y",
    "TCA" => "S",
    "TCC" => "S",
    "TCG" => "S",
    "TCT" => "S",
    "TGA" => "*",
    "TGC" => "C",
    "TGG" => "W",
    "TGT" => "C",
    "TTA" => "L",
    "TTC" => "F",
    "TTG" => "L",
    "TTT" => "F",
    }
}

fn dicodons_by_amino_acids() -> HashMap<String, Vec<String>> {
    let mut res = HashMap::new();
    for (codon1, aa1) in codon_to_aa_map().iter() {
        for (codon2, aa2) in codon_to_aa_map().iter() {
            let dicodon = format!("{}{}", codon1, codon2);
            let key = format!("{}{}", aa1, aa2);
            res.entry(key).or_insert(Vec::new()).push(dicodon);
        }
    }
    res
}

#[derive(Serialize)]
pub struct AnalyzeResult {
    pub scores: Vec<f64>,
    pub codons: Vec<String>,
    pub amino_acids: Vec<String>,
    pub total_score: f64
}

impl AnalyzeResult {
    fn new(
        dna_sequence: &str,
        dicodon_frequencies: &HashMap<String, f64>,
    ) -> Result<AnalyzeResult, String> {
        if dna_sequence.len() % 3 != 0 {
            return Err("DNA length was not a multiple of three - must be full codons".into());
        }
        let mut scores: Vec<f64> = Vec::new();
        let mut codons: Vec<String> = Vec::new();
        let mut amino_acids: Vec<String> = Vec::new();
        let mut last_codon = "";
        let mut total_score = 0.0;
        let codon_map = codon_to_aa_map();
        for (first_codon, second_codon) in dna_sequence.as_bytes().chunks(3).tuple_windows() {
            let first_codon = std::str::from_utf8(first_codon).unwrap();
            let second_codon = std::str::from_utf8(second_codon).unwrap();
            let dicodon = format!("{}{}", first_codon, second_codon);
            let f = dicodon_frequencies.get(&dicodon).unwrap_or(&0.);
            scores.push(*f);
            total_score += *f;
            codons.push(first_codon.to_string());
            amino_acids.push(codon_map.get(first_codon).unwrap().to_string());
            last_codon = second_codon;
        }
        codons.push(last_codon.to_string());
        let last_aa: String = codon_map.get(last_codon).unwrap().to_string();
        amino_acids.push(last_aa);
        Ok(AnalyzeResult {
            scores,
            codons,
            amino_acids,
            total_score,
        })
    }
}

#[wasm_bindgen]
pub fn analyze(dicodon_frequency_table: &str, sequence: &str) -> Result<JsValue, JsValue> {
    let dicodon_frequencies = parse_dicodon_frequency_table(dicodon_frequency_table)?;
    let parsed = parse_sequence(sequence)?;
    match parsed {
        SequenceKind::DNA(s) => {
            Ok(JsValue::from_serde(&AnalyzeResult::new(&s, &dicodon_frequencies).unwrap()).unwrap())
        }
        SequenceKind::AA(_) => Err(
            "Can only analyze DNA codon, AA sequences can only be optimized \
            (any N is autodetected as Asparagin)"
                .into(),
        ),
    }
}

#[derive(Serialize)]
pub struct OptimizeResult {
    pub optimized_sequence: String,
    pub analyzed: AnalyzeResult,
    pub analyzed_input: Option<AnalyzeResult>,
}

#[wasm_bindgen]
pub fn optimize(dicodon_frequency_table: &str, sequence: &str) -> Result<JsValue, JsValue> {
    //these are AGTCCT -> probability
    let dicodon_frequencies = parse_dicodon_frequency_table(dicodon_frequency_table)?;
    let parsed = parse_sequence(sequence)?;
    log!("parsed");
    let (aa_seq, analyzed_input) = match parsed {
        SequenceKind::DNA(s) => {
            if s.len() % 3 != 0 {
                return Err("DNA length was not a multiple of three - must be full codons".into());
            }
            let map = codon_to_aa_map();
            let mut aa = "".to_string();
            for ii in (0..s.len()).step_by(3) {
                log!("{}, {}", ii, s.len());
                let codon = &s[ii..ii + 3];
                aa += map.get(codon).expect("Invalid codon found?");
            }
            (aa, Some(AnalyzeResult::new(&s, &dicodon_frequencies).unwrap()))
        }
        SequenceKind::AA(s) => (s,None),
    };
    let states = codon_to_aa_map()
        .into_iter()
        .map(|(k, v)| (k, v.as_bytes()[0]))
        .collect();
    //these rae ("AGT","CCT" -> probability
    let transitions = dicodon_frequencies
        .iter()
        .map(|(k, v)| ((k[0..3].to_string(), k[3..6].to_string()), *v))
        .collect();
    let hmm = DegenerateHMM::new(states, transitions).expect("Failed to build hmm");
    let (_max_score, max_sequence) = hmm.viterbi(aa_seq.as_bytes()).expect("Could not find path");
    let optimized_sequence = max_sequence.join("");
    let analyzed = AnalyzeResult::new(&optimized_sequence, &dicodon_frequencies).unwrap();
    let res = OptimizeResult {
        optimized_sequence,
        analyzed,
        analyzed_input,
    };
    Ok(JsValue::from_serde(&res).unwrap())
}

fn parse_dicodon_frequency_table(input: &str) -> Result<HashMap<String, f64>, String> {
    let mut res = HashMap::new();
    let lines = input.trim().split("\n");
    let mut lines = lines.filter(|x| !x.starts_with("#"));
    let header = lines
        .next()
        .ok_or("Dicodonfrequencies contained no lines")?;
    if header != "AA\tDicodon\tRel.Freq\tCount" {
        return Err(format!(
            "Header was not 'AADicodon Rel.Freq Count' was '{:?}'",
            header
        ));
    }
    for (ii, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let mut parts = line.split("\t");
        let _aa = parts
            .next()
            .ok_or_else(|| format!("frequency table: error handling line {}, failed to find AA", ii+1))?;
        let dicodon = parts
            .next()
            .ok_or_else(|| format!("frequency table: error handling line {}, failed to find dicodon", ii+1))?
            .to_uppercase();
        let rel_freq = parts
            .next()
            .ok_or_else(|| format!("frequency table: error handling line {}, failed to find rel_freq", ii+1))?;
        let _count = parts
            .next()
            .ok_or_else(|| format!("frequency table: error handling line {}, failed to find count", ii+1))?;
        let rel_freq: f64 = rel_freq.parse::<f64>().or_else(|_| {
            Err(format!(
                "frequency table: error handling line {}, failed to parse relative frequency as a float",
                ii+1
            ))
        })?;
        res.insert(dicodon, rel_freq);
    }
    Ok(res)
}

enum SequenceKind {
    DNA(String),
    AA(String),
}

fn parse_sequence(input: &str) -> Result<SequenceKind, String> {
    let input = input.trim().to_uppercase();
    let lines = input
        .split("\n")
        .skip(if input.starts_with(">") { 1 } else { 0 });
    let mut sequence = "".to_string();
    let mut all_dna = true;
    for l in lines {
        for c in l.chars() {
            match c {
                'A' | 'C' | 'G' | 'T' => sequence.push(c),
                //'A' | 'C' | '| 'G'  'T'   are also valid aa, but handled bove
                'R' | 'N' | 'D' | 'B' | 'E' | 'Q' | 'Z' | 'H' | 'I' | 'L' | 'K' | 'M' | 'F'
                | 'P' | 'S' | 'W' | 'Y' | 'V' => {
                    all_dna = false;
                    sequence.push(c)
                }
                ' ' | '\n' | '\t' => {}
                '>' => {
                    return Err(format!(
                        "Invalid letter in sequence: '{}'. Must be single fasta entry",
                        c
                    ))
                }
                _ => return Err(format!("Invalid letter in sequence: '{}'", c)),
            };
        }
    }
    Ok(match all_dna {
        true => SequenceKind::DNA(sequence),
        false => SequenceKind::AA(sequence),
    })
}

/// A HMM that only has one emission on each state
struct DegenerateHMM {
    states: HashMap<String, u8>,
    transition_scores: HashMap<(String, String), f64>,
    allowed_emissions: HashSet<u8>,
}

impl DegenerateHMM {
    fn new(
        states: HashMap<String, u8>,
        transition_probabilities: HashMap<(String, String), f64>,
    ) -> Result<DegenerateHMM, String> {
        let mut transition_scores = HashMap::new();
        for (k, v) in transition_probabilities.into_iter() {
            if (v < 0.) || (v > 1.) {
                return Err("Non probability in transition probabilities".to_string());
            }
            if !states.contains_key(&k.0) {
                return Err("Invalid from state".to_string());
            }
            if !states.contains_key(&k.1) {
                return Err("Invalid to state".to_string());
            }

            transition_scores.insert(k, v.ln());
        }
        let allowed_emissions = states.values().map(|x| *x).collect();
        Ok(DegenerateHMM {
            states,
            transition_scores,
            allowed_emissions,
        })
    }

    /// Given a list of observed emissions,
    /// give us the maximal score and path through the hidden markov model
    /// @observed_sequence: an iterable of symbols
    /// result: (maximum_score, [state_0, state_1, state_2..., state_len(observed_sequence)])

    fn viterbi(&self, observed_sequence: &[u8]) -> Result<(f64, Vec<String>), String> {
        if observed_sequence.is_empty() {
            return Ok((0., Vec::new()));
        }
        #[allow(non_snake_case)]
        let mut V = Vec::new(); // this could be much faster if we mapped the states to usize and used an array instead.
        for _ in 0..observed_sequence.len() {
            V.push(HashMap::new())
        }
        let mut path = HashMap::new();
        for (state, emission) in self.states.iter() {
            V[0].insert(
                state,
                if *emission == observed_sequence[0] {
                    0.
                } else {
                    f64::NEG_INFINITY
                },
            );
            //V[0][state] = 0 if emission == observed_sequence[0] else float("-inf")
            path.insert(state, vec![state]);
        }

        for t in 1..observed_sequence.len() {
            let mut new_path = HashMap::new(); // we only keep the survivor paths for the latest t
            let mut any_survivors = false;
            for (to_state, emission) in self.states.iter() {
                // # for each state we could potentially end up
                if *emission == observed_sequence[t] {
                    //   this is valid transition, we can end up here
                    // let's find out what's the most likely state we could have come from
                    let mut survivor_score = f64::NEG_INFINITY;
                    let mut survivor_state = None;
                    for from_state in self.states.keys() {
                        //how expensive is it to go from from_state to to_state?
                        let transition_score = self
                            .transition_scores
                            .get(&(from_state.to_string(), to_state.to_string()))
                            .expect("invalid transition");
                        // and how good was being in from_state in the first place
                        let score_for_being_in_from_state = V[t - 1].get(from_state).unwrap();
                        let score_here = transition_score + score_for_being_in_from_state; // + in log is * in the normal scale
                        if score_here > survivor_score {
                            // it's best (so far) to have come from from_state!
                            survivor_score = score_here;
                            survivor_state = Some(from_state);
                            any_survivors = true; // this Hmm can have produced this sequence
                        }
                    }
                    match survivor_state {
                        Some(survivor_state) => {
                            let mut ns = path[survivor_state].clone();
                            ns.push(to_state);
                            new_path.insert(to_state, ns);
                            V[t].insert(to_state, survivor_score);
                        }
                        None => {
                            V[t].insert(to_state, f64::NEG_INFINITY);
                        }
                    }
                } else {
                    V[t].insert(to_state, f64::NEG_INFINITY);
                }
            }

            if !any_survivors {
                return Err(format!(
                    "This HMM can't have generated this sequence in position {}",
                    t
                ));
            }
            path = new_path
        }
        let final_scores = &V[V.len() - 1]; // hashmap state->score where we could have ended up
        let mut maximum_final_score = std::f64::NEG_INFINITY;
        let mut maximum_final_state = None;
        for state in final_scores.keys() {
            let fs = *final_scores.get(state).unwrap();
            if fs > maximum_final_score {
                maximum_final_score = fs;
                maximum_final_state = Some(state);
            }
        }
        match maximum_final_state {
            None => Err(format!(
                "This HMM can't have generated this sequence in position {}",
                observed_sequence.len() - 1
            )),
            Some(mf) => Ok((
                std::f64::consts::E.powf(maximum_final_score),
                path.get(mf)
                    .unwrap()
                    .into_iter()
                    .map(|x| x.to_string())
                    .collect(),
            )),
        }

        /*
                    """
                if not observed_sequence:
                    return 0, []

                V = (
                    []
                )  # the score of the survivor path that ended in state y at observed_sequence t
                for t in range(0, len(observed_sequence)):
                    if observed_sequence[t] not in self.allowed_emissions:
                        raise ValueError("Invalid emission: %s" % observed_sequence[t])
                    V.append({})

                path = (
                    {}
                )  # for the position we have reached in observed_sequence, the survivor paths for each possible state

                # initialize for leaving the 'virtual' beginning state'
                # score is 0 for a valid transition,
                # so this 'virtual' begin state does not affect the final score
                for state, emission in self.states.items():
                    V[0][state] = 0 if emission == observed_sequence[0] else float("-inf")
                    path[state] = [state]

                for t in range(1, len(observed_sequence)):
                    new_path = {}  # we only keep the survivor paths for the latest t

                    any_survivors = (
                        False  # to check if this HMM could have reached this state...
                    )
                    for (
                        to_state,
                        emission,
                    ) in self.states.items():  # for each state we could potentially end up
                        if (
                            emission == observed_sequence[t]
                        ):  # this is valid transition, we can end up here
                            # let's find out what's the most likely state we could have come from
                            survivor_score = float("-inf")
                            survivor_state = None
                            for from_state, emission in self.states.items():
                                # how expensive is it to go from from_state to to_state?
                                transition_score = self.transition_scores[from_state, to_state]
                                # and how good was being in from_state in the first place
                                score_for_being_in_from_state = V[t - 1][from_state]
                                score_here = (
                                    transition_score + score_for_being_in_from_state
                                )  # + in log is * in the normal scale
                                if (
                                    score_here > survivor_score
                                ):  # it's best (so far) to have come from from_state!
                                    survivor_score = score_here
                                    survivor_state = from_state
                                    any_survivors = (
                                        True  # this Hmm can have produced this sequence
                                    )
                            if (
                                survivor_state is not None
                            ):  # there is no way we can have ended up in to_state. But there might be other states with the right emission
                                new_path[to_state] = path[survivor_state] + [to_state]
                                V[t][to_state] = survivor_score
                            else:
                                V[t][to_state] = float("-inf")
                        else:
                            # new_path[to_state] = None #no access to this path should happen. Not having this line would make it throw a KeyError if access did happen.
                            V[t][to_state] = float("-inf")

                    if not any_survivors:
                        raise ModelSequenceMismatchException(
                            "This HMM can't have generated this sequence in position %i" % t
                        )
                    path = new_path
                final_scores = V[-1]
                maximum_final_score = float("-inf")
                maximum_final_state = None
                for state in final_scores:
                    if final_scores[state] > maximum_final_score:
                        maximum_final_score = final_scores[state]
                        maximum_final_state = state
                if maximum_final_state is None:
                    raise ModelSequenceMismatchException(
                        "This HMM can't have generated this sequence in position %i"
                        % len(observed_sequence)
                    )
                return pow(2, maximum_final_score), path[maximum_final_state]


        &*/
    }

    fn score(&self, observed_sequence: &[&str]) -> Result<f64, String> {
        if observed_sequence.len() == 1 {
            return Ok(0.0);
        }
        let mut score = 0.0;
        for t in 1..observed_sequence.len() {
            let key = (
                observed_sequence[t - 1].to_string(),
                observed_sequence[t].to_string(),
            );
            match self.transition_scores.get(&key) {
                Some(x) => score += x,
                None => {
                    return Err("Model can not have generated this sequence".to_string());
                }
            }
        }
        Ok(std::f64::consts::E.powf(score))
    }
}

#[cfg(test)]
mod test {
    use super::DegenerateHMM;
    use std::collections::HashMap;

    #[test]
    fn test_ok() {
        let states = collection! {
        "a" => b'A',
        "b" => b'B',};
        let transitions = collection! {
            ("a".to_string(),"a".to_string()) => 0.6,
            ("b".to_string(),"a".to_string()) => 0.5,
            ("a".to_string(),"b".to_string()) => 0.5,
        };

        assert!(DegenerateHMM::new(states, transitions).is_ok());
    }
    #[test]
    fn test_non_prop() {
        let states: HashMap<String, u8> = collection! {
        "a" => b'A',
        "b" => b'B',};
        let transitions = collection! {
            ("a".to_string(),"a".to_string()) => 0.6,
            ("b".to_string(),"a".to_string()) => 0.5,
            ("a".to_string(),"b".to_string()) => 1.1,
        };
        assert!(DegenerateHMM::new(states.clone(), transitions).is_err());
        let transitions = collection! {
            ("a".to_string(),"a".to_string()) => 0.6,
            ("b".to_string(),"a".to_string()) => -0.5,
            ("a".to_string(),"b".to_string()) => 0.1,
        };
        assert!(DegenerateHMM::new(states, transitions).is_err());
    }

    use approx::relative_eq;

    #[test]
    fn test_scores() {
        let states: HashMap<String, u8> = collection! {
        "a" => b'A',
        "b" => b'B',};
        let transitions = collection! {
            ("a".to_string(),"b".to_string()) => 0.4,
            ("a".to_string(),"a".to_string()) => 0.6,
            ("b".to_string(),"a".to_string()) => 0.5,
            ("b".to_string(),"b".to_string()) => 0.5,
        };
        let hmm = DegenerateHMM::new(states.clone(), transitions).unwrap();
        assert!(hmm.score(&vec!["a"][..]).unwrap() == 0.0); //no score for entering the sequence
        assert!(hmm.score(&vec!["b"][..]).unwrap() == 0.0);
        assert!(hmm.score(&vec!["a", "b"][..]).unwrap() == 0.4);
        assert!(hmm.score(&vec!["b", "a"][..]).unwrap() == 0.5);
        println!("was {}", hmm.score(&vec!["a", "b", "b"][..]).unwrap());

        assert!(relative_eq!(
            hmm.score(&vec!["a", "b", "b"][..]).unwrap(),
            0.4 * 0.5
        ));
        assert!(relative_eq!(
            hmm.score(&vec!["b", "a", "b"][..]).unwrap(),
            0.5 * 0.4
        ));
        assert!(relative_eq!(
            hmm.score(&vec!["b", "a", "a"][..]).unwrap(),
            0.5 * 0.6
        ));
    }

    #[test]
    fn test_init_raises_on_invalid_from_state() {
        let states: HashMap<String, u8> = collection! {
        "a" => b'A',
        "b" => b'B',};
        let transitions = collection! {
            ("X".to_string(),"b".to_string()) => 0.4,
            ("a".to_string(),"a".to_string()) => 0.6,
            ("b".to_string(),"a".to_string()) => 0.5,
            ("b".to_string(),"b".to_string()) => 0.5,
        };
        assert!(DegenerateHMM::new(states.clone(), transitions).is_err());
    }
    #[test]
    fn test_init_raises_on_invalid_to_state() {
        let states: HashMap<String, u8> = collection! {
        "a" => b'A',
        "b" => b'B',};
        let transitions = collection! {
            ("a".to_string(),"b".to_string()) => 0.4,
            ("a".to_string(),"a".to_string()) => 0.6,
            ("b".to_string(),"a".to_string()) => 0.5,
            ("b".to_string(),"X".to_string()) => 0.5,
        };
        assert!(DegenerateHMM::new(states.clone(), transitions).is_err());
    }
    fn get_trivial_hmm() -> DegenerateHMM {
        let states = collection! {"raining"=> b'h', // stay at home
                "sunny"=> b'w' // go for a walk
        };
        let transition_scores = collection! {
        ("raining".to_string(), "raining".to_string())=> 0.7,
        ("raining".to_string(), "sunny".to_string())=> 0.3,
        ("sunny".to_string(), "raining".to_string())=> 0.9,
        ("sunny".to_string(), "sunny".to_string())=> 0.1,
        }; //yeah... it rains a lot
        DegenerateHMM::new(states, transition_scores).unwrap()
    }
    #[test]
    fn test_trivial_viterbi() {
        let hmm = get_trivial_hmm();
        let observed = &[b'h', b'h', b'w', b'w', b'h'];
        let (score, state_sequence) = hmm.viterbi(observed).unwrap();
        let supposed = vec!["raining", "raining", "sunny", "sunny", "raining"];
        assert!(state_sequence == supposed);
        assert!(relative_eq!(score, 0.7 * 0.3 * 0.1 * 0.9));
        assert!(relative_eq!(score, hmm.score(&supposed[..]).unwrap()));
    }

    #[test]
    fn test_raises_on_invalid_emission_in_viterbi() {
        let hmm = get_trivial_hmm();
        let observed = &[b'h', b'x', b'w', b'w', b'h'];
        assert!(hmm.viterbi(observed).is_err());
    }

    #[test]
    fn test_raises_on_imposible_sequence_to_score() {
        let states = collection! {"raining"=> b'h',
        "sunny"=> b'w'};
        let transition_scores = collection! {
        ("raining".to_string(), "raining".to_string())=> 0.7,
        ("sunny".to_string(), "raining".to_string())=> 0.9,
        ("sunny".to_string(), "sunny".to_string())=> 0.1,
        }; //yeah... it rains a lot
        let hmm = DegenerateHMM::new(states, transition_scores).unwrap();
        assert!(hmm.score(&["sunny", "raining", "sunny"]).is_err()) //# no transition raining->sunny
    }

    #[test]
    fn test_raises_on_invalid_state() {
        let states = collection! {"raining"=> b'h',
        "sunny"=> b'w'};
        let transition_scores = collection! {
        ("raining".to_string(), "raining".to_string())=> 0.7,
        ("raining".to_string(), "sunny".to_string())=> 0.7,
        ("sunny".to_string(), "raining".to_string())=> 0.9,
        ("sunny".to_string(), "sunny".to_string())=> 0.1,
        }; //yeah... it rains a lot
        let hmm = DegenerateHMM::new(states, transition_scores).unwrap();
        assert!(hmm.score(&["sunny", "snowstorm", "sunny"]).is_err()) //# no snowstorm state
    }

    #[test]
    fn test_empty_sequence_viterbi() {
        let hmm = get_trivial_hmm();
        let (max_score, state_sequence) = hmm.viterbi(&[]).unwrap();
        assert!(max_score == 0.0);
        assert!(state_sequence.len() == 0);
    }

    #[test]
    fn test_multiple_states_same_emission() {
        let states = collection! {"raining"=> b'h',
        "sunny"=> b'w',
        "snowstorm"=> b'h'};
        let transition_scores = collection! {
        ("raining".to_string(), "raining".to_string())=> 0.7,
        ("raining".to_string(), "sunny".to_string())=> 0.3,
        ("sunny".to_string(), "raining".to_string())=> 0.9,
        ("sunny".to_string(), "sunny".to_string())=> 0.1,
        ("snowstorm".to_string(), "sunny".to_string())=> 0.4,
        ("snowstorm".to_string(), "snowstorm".to_string())=> 0.5,
        ("snowstorm".to_string(), "raining".to_string())=> 0.1,
        ("raining".to_string(), "snowstorm".to_string())=> 0.1,
        ("sunny".to_string(), "snowstorm".to_string())=> 0.2};
        let hmm = DegenerateHMM::new(states, transition_scores).unwrap();
        let observed = &[b'w', b'h', b'h', b'h', b'w'];
        let (max_score, state_sequence) = hmm.viterbi(observed).unwrap();
        let supposed = &["sunny", "raining", "raining", "raining", "sunny"];
        let supposed_score = hmm.score(supposed).unwrap();
        assert!(state_sequence == supposed);
        assert!(relative_eq!(max_score, supposed_score));
    }
    #[test]
    fn test_multiple_states_same_emission_make_it_snow() {
        let states = collection! {"raining"=> b'h',
        "sunny"=> b'w',
        "snowstorm"=> b'h'
        };
        let transition_scores = collection! {
        ("raining".to_string(), "raining".to_string())=> 0.1,
        ("raining".to_string(), "sunny".to_string())=> 0.8,
        ("raining".to_string(), "snowstorm".to_string())=> 0.1,
        ("sunny".to_string(), "raining".to_string())=> 0.4,
        ("sunny".to_string(), "sunny".to_string())=> 0.5,
        ("sunny".to_string(), "snowstorm".to_string())=> 0.2, //snow storm"s seldom start
        ("snowstorm".to_string(), "sunny".to_string())=> 0.100000001, //otherwise it's a tie <F12
        ("snowstorm".to_string(), "snowstorm".to_string())=> 0.8, //#but they keep on for a long time
        ("snowstorm".to_string(), "raining".to_string())=> 0.1};
        let hmm = DegenerateHMM::new(states, transition_scores).unwrap();
        let observed = [b'w', b'h', b'w'];
        let (max_score, state_sequence) = hmm.viterbi(&observed).unwrap();
        let supposed = ["sunny", "raining", "sunny"]; //a single day is most likely just raining
        let supposed_score = hmm.score(&supposed).unwrap();
        assert!(relative_eq!(max_score, supposed_score));

        let observed = [b'w', b'h', b'h', b'w'];
        let (max_score, state_sequence) = hmm.viterbi(&observed).unwrap();
        let supposed = ["sunny", "raining", "raining", "sunny"]; //#two days: still raining
        let supposed_score = hmm.score(&supposed).unwrap();
        assert!(state_sequence == supposed);
        assert!(relative_eq!(max_score, supposed_score));

        let observed = [b'w', b'h', b'h', b'h', b'w'];
        let (max_score, state_sequence) = hmm.viterbi(&observed).unwrap();
        let supposed = ["sunny", "snowstorm", "snowstorm", "snowstorm", "sunny"]; //#while three days in a row is snowing
        let supposed_score = hmm.score(&supposed).unwrap();
        println!("{:?}", state_sequence);
        assert!(state_sequence == supposed);
        assert!(relative_eq!(max_score, supposed_score));

        let observed = [b'w', b'h', b'h', b'h', b'w', b'h'];
        let (max_score, state_sequence) = hmm.viterbi(&observed).unwrap();
        let supposed = [
            "sunny",
            "snowstorm",
            "snowstorm",
            "snowstorm",
            "sunny",
            "raining",
        ]; //and some more rain
        let supposed_score = hmm.score(&supposed).unwrap();
        assert!(state_sequence == supposed);
        assert!(relative_eq!(max_score, supposed_score));
    }
}
