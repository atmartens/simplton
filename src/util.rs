/***
SIMP'LTON: simulation of polysome translation
    Copyright (C) 2025 Andrew T. Martens

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License, version 3,
    as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

e-mail: simplton@logical.bio
***/

use anyhow::{Context, Result, bail};
use std::{fs::File, io::{BufRead, BufReader}, path::{Path, PathBuf}};
use csv::ReaderBuilder;
use std::cmp::Ordering;
use bio::io::fasta::{Reader, Record};

/// Returns `true` if the given sequence (ASCII‐encoded) is composed *exclusively*
/// of the letters A, T, G, C, or U (case‐insensitive).  
/// In other words, it returns `true` if NO character falls outside [AaTtGgCcUu].
fn looks_like_nucleotide_only(seq: &[u8]) -> bool {
    seq
	.iter()
	.all(
	    |&b| matches!(
		b,
		b'A' |
		b'a' |
		b'T' |
		b't' |
		b'G' |
		b'g' |
		b'C' |
		b'c' |
		b'U' |
		b'u'
	    )
	)
}

pub fn read_first_protein_sequence(path : &str) -> Result<String> {
    let file = File::open(path)
	.with_context(|| format!("Failed to open FASTA file '{}'", path))?;
    let reader = BufReader::new(file);
    let fasta_reader = Reader::new(reader);

    if let Some(next) = fasta_reader.records().next() {
	let record: Record = next.expect("Error while reading the first FASTA record");
	let id = record.id();
	let seq_bytes = record.seq();
	
	if looks_like_nucleotide_only(seq_bytes) {
	    eprintln!(
		"Warning: sequence {} contains only nucleotides. \
		 Are you sure it's DNA and not protein?", id
	    )
	}
	
	let seq_string = String::from_utf8(seq_bytes.to_vec())
	    .context("Invalid UTF-8 found in sequence data")?;
	
	Ok(seq_string)
    } else {
	bail!("FASTA file contained no records")
    }

}

/* Input a file containing rates for each codon. The file is formatted with tab-
separated fields, with the first element denoting the position, and the second
the rate constant. One pair per line. */
pub fn init_rates_from_file(filename: &str) -> Result<Vec<f32>> {
    let file = BufReader::new(File::open(PathBuf::from(&filename))?);
    let mut p = Vec::new();
    // NOTE: start codon is at zero!
    for l in file.lines() {
        let l = l?;
        if !l.starts_with("#") {
            let a: Vec<_> = l.split("\t").collect();
            let position = a[0].parse()?;
	    
            let rate = a[1].parse()?;
	    if rate < 0.0 {
		bail!("Invalid translation rate detected at position {}", position);
	    }
            p.insert(position, rate);
        }
    }

    Ok(p)
}

/* Input a string describing the rates at each codon, return a vector of rate 
constants. */
pub fn parse_rates(input: &str) -> Result<Vec<f32>> {
    let mut r = Vec::new();
    let segments: Vec<_> = input.split(",").collect();

    for s in segments {
        let e: Vec<_> = s.split(":").collect();
        match e.len() {
            3 => {
                let rate  = e[0].parse()?;
                let begin = e[1].parse::<usize>()?;
                let end   = e[2].parse::<usize>()?;

                if (rate < 0.0) || (begin < 1) || (end < 1) || (begin > end) {
                    bail!("conflicting rate information in rates string!");
                } else {
                    for _ in begin..end + 1 {
                        r.push(rate);
                    }
                }
            }

            2 => {
                let rate  = e[0].parse()?;
                let begin = e[1].parse::<usize>()?;

                if (rate < 0.0) || (begin < 1) {
                    bail!("conflicting rate information in rates string!");
                } else {
                    r.push(rate);
                }
            }

            _ => bail!("non-numeric values encountered in rates string."),
        }
    }

    Ok(r)
}

pub fn init_events_tsv<P: AsRef<Path>>(
    path: P
) -> Result<Vec<(f32, Vec<usize>, Vec<f32>)>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_path(path)?;

    let mut events = Vec::new();
    for result in rdr.records() {
        let record = result?;
        let time: f32 = record[0].parse()?;

	let positions = record[1]
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<usize>())
            .collect::<Result<Vec<_>, _>>()?;

	let rates = record[2]
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<f32>())
            .collect::<Result<Vec<_>, _>>()?;
	
	if positions.len() != rates.len() {
	    bail!(
		"mismatched lengths at time {}: {} positions vs {} rates",
		time, positions.len(), rates.len()
	    );
}

	events.push((time, positions, rates));
    }

    events.sort_by(|a,b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));
    
    Ok(events)
}


/// Parse a semicolon-separated list of events of the form
///    time:pos1,pos2,…:rate1,rate2,…
/// into a Vec of (time, positions, rates), sorted descending by time
/// so that `pop()` yields the next‐earliest event.
///
/// e.g.    
/// ```ignore
/// let s = "1.0:1,2,3:0.1,0.2,0.3;\
///          0.5:4:0.4;\
///          2.0:5,6:0.5,0.6;";
/// let mut ev = parse_events_str(s)?;
/// assert_eq!(ev.pop().unwrap().0, 0.5); // 0.5 is the smallest time
/// ```
pub fn parse_events_str(input: &str) -> Result<Vec<(f32, Vec<usize>, Vec<f32>)>> {
    let mut events = Vec::new();

    for (idx, chunk) in input.split(';').enumerate() {
        let chunk = chunk.trim();
        if chunk.is_empty() {
            continue;
        }

        let parts: Vec<&str> = chunk.split(':').collect();
        if parts.len() != 3 {
            bail!(
                "expected 3 ':'-separated fields in event {}: `{}`",
                idx + 1,
                chunk
            );
        }

        // parse the timestamp
        let time: f32 = parts[0].parse()
            .with_context(|| format!("invalid time in event {}: `{}`", idx + 1, parts[0]))?;

        // parse positions
        let positions = parts[1]
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<usize>())
            .collect::<Result<Vec<_>, _>>()?;

        // parse rates
        let rates = parts[2]
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<f32>())
            .collect::<Result<Vec<_>, _>>()?;

        // ensure they line up
        if positions.len() != rates.len() {
            bail!(
                "length mismatch in event {}: {} positions vs {} rates",
                idx + 1,
                positions.len(),
                rates.len()
            );
        }

        events.push((time, positions, rates));
    }

    // sort descending by time -> pop() gives smallest/earliest
    events.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Ordering::Equal));

    Ok(events)
}

