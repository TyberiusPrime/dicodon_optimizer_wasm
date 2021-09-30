#![allow(dead_code)]
mod utils;

//use seq_io::fasta::{Reader, Record};
use itertools::Itertools;
use std::collections::HashMap;
use std::convert::TryInto;
use wasm_bindgen::prelude::*;
use serde::{Serialize, Deserialize};

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
}

#[wasm_bindgen]
pub fn analyze(dicodon_frequency_table: &str, sequence: &str) -> Result<JsValue, JsValue> {
    let dicodon_frequencies = parse_dicodon_frequency_table(dicodon_frequency_table)?;
    let parsed = parse_sequence(sequence)?;
    match parsed {
        SequenceKind::DNA(s) => {
            if s.len() % 3 != 0 {
                return Err("DNA length was not a multiple of three - must be full codons".into());
            }
            let mut scores: Vec<f64> = Vec::new();
            let mut codons: Vec<String> = Vec::new();
            let mut amino_acids: Vec<String> = Vec::new();
            let mut last_codon = "";
            let codon_map = codon_to_aa_map();
            for (first_codon, second_codon) in s.as_bytes().chunks(3).tuple_windows() {
                let first_codon = std::str::from_utf8(first_codon).unwrap();
                let second_codon = std::str::from_utf8(second_codon).unwrap();
                let dicodon = format!("{}{}", first_codon, second_codon);
                let f = dicodon_frequencies.get(&dicodon).unwrap_or(&0.);
                scores.push(*f);
                codons.push(first_codon.to_string());
                amino_acids.push(codon_map.get(first_codon).unwrap().to_string());
                last_codon = second_codon;
            }
            codons.push(last_codon.to_string());
            let last_aa: String = codon_map.get(last_codon).unwrap().to_string();
            amino_acids.push(last_aa);
            Ok(JsValue::from_serde(&AnalyzeResult {
                scores,
                codons, 
                amino_acids,
            }).unwrap())
        }
        SequenceKind::AA(s) => {
            Err("Can only analyze DNA codon, AA sequences can only be optimized".into())
        }
    }
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
        let mut parts = line.split("\t");
        let _aa = parts
            .next()
            .ok_or_else(|| format!("error handling line {}, failed to find AA", ii))?;
        let dicodon = parts
            .next()
            .ok_or_else(|| format!("error handling line {}, failed to find dicodon", ii))?
            .to_uppercase();
        let rel_freq = parts
            .next()
            .ok_or_else(|| format!("error handling line {}, failed to find rel_freq", ii))?;
        let _count = parts
            .next()
            .ok_or_else(|| format!("error handling line {}, failed to find count", ii))?;
        let rel_freq: f64 = rel_freq.parse::<f64>().or_else(|_| {
            Err(format!(
                "error handling line {}, failed to parse relative frequency as a float",
                ii
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
