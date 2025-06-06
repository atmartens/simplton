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

/* main.rs: command-line interface */

#![allow(non_snake_case)] // allow variables with "mRNA" in them

use polars::prelude::{ParquetReader, ParquetWriter, SerReader};

use threadpool::ThreadPool;

use std::fs::{create_dir_all, File, write, remove_file};
use std::io::{self, BufWriter, Write};
use std::path::{PathBuf, Path};

use rand_distr::{Distribution, Exp};

use clap::{Parser, Subcommand, ArgAction};

use arrow2::datatypes::{Schema, Field, DataType};
use arrow2::array::*;
use arrow2::chunk::Chunk;
use arrow2::io::parquet::write::{
    CompressionOptions, WriteOptions, FileWriter, Version, Encoding,
    RowGroupIterator
};

use anyhow::{Context, Result, bail};

use std::sync::mpsc::channel;
use std::sync::Arc;
use std::env;

mod simplton;
mod dataframe_operations;
mod util;

/* Simple struct for receiving results from a thread */
pub struct SimulationMessage {
    mRNA_number: u16,
    trajectories: Result<Chunk<Arc<dyn Array>>>,
}

#[derive(Parser)]
#[command(
    name = "simplton",
    version = "1.0.0",
    author = "Andrew T. Martens <andrew_martens@hms.harvard.edu>",
    about = "SIMP'LTON -- SIMulating PoLysome TranslatiON",
    disable_help_subcommand = true,
    disable_version_flag = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(
	about = "Simulate ribosome movement and collisions, given the \
		 translation rates at each codon"
    )]
    Simulate {
        #[arg(
	    short = 't',
	    long = "simulation-time",
	    help = "Time in seconds to run the simulation",
	    required = true
	)]
        simulation_time: f32,

        #[arg(
	    short = 'o',
	    long = "out-dir",
	    help = "Directory to write to",
	    required = false
	)]
        out_dir: Option<String>,

        #[arg(
	    short = 'λ',
	    long = "mRNA-half-life",
	    default_value_t = f32::INFINITY,
	    help = "mRNA half-life, in seconds",
	    required = false
	)]
        mRNA_half_life: f32,

        #[arg(
	    short = 'n',
	    long = "num-mRNA",
	    required_unless_present = "mRNA_freq",
	    conflicts_with = "mRNA_freq",
	    help = "If fixed mRNA, how many?"
	)]
        num_mRNA: Option<u16>,

        #[arg(
	    short = 'τ',
	    long = "mRNA-freq",
	    required_unless_present = "num_mRNA",
	    conflicts_with = "num_mRNA",
	    help = "If spawned mRNA, how often?",
	)]
        mRNA_freq: Option<f32>,

	// Rates
        #[arg(
	    short = 'r',
	    long = "rates-string",
	    required_unless_present = "rates_file",
	    conflicts_with = "rates_file",
	    help = "Parse the translation rates from the command line"
	)]
        rates_string: Option<String>,

        #[arg(
	    short = 'f',
	    long = "rates-file",
	    required_unless_present = "rates_string",
	    conflicts_with = "rates_string",
	    help = "Input the translation rates from a file",
	)]
        rates_file: Option<String>,

	// Events
	#[arg(
	    short = 'E',
	    long = "events-string",
	    help = "Input the events from a string",
	    conflicts_with = "events_file"
	)]
        events_string: Option<String>,
	
	#[arg(
	    short = 'e',
	    long = "events-file",
	    help = "Input the events from a file",
	    conflicts_with = "events_string"
	)]
        events_file: Option<String>,

        #[arg(
	    short = 'T',
	    long = "num-threads",
	    default_value_t = 1,
	    help = "Maximum number of threads"
	)]
        num_threads: usize,

	#[arg(short, long, action = ArgAction::Count)]
	verbosity: u8,

	#[arg(
	    short = 'L',
	    long = "footprint-length",
	    default_value_t = 10,
	    help = "Size of a ribosome footprint, in codons"
	)]
        footprint_length: usize,

	#[arg(
	    short = 'A',
	    long = "a-pos-length",
	    default_value_t = 4,
	    help = "Location of A-site position, in codons"
	)]
        a_pos_length: usize,
    },

    #[command(
	about = "Calculate statistics of a translation simulation: step times,\
		 final occupancy, protein production (both by mRNA and total), \
		 synthesis times and collision frequencies."
    )]
    Statistics {
        #[arg(
	    short = 't',
	    long = "trajectories_path",
	    default_value = "trajectories.pq",
	    help = "Trajectories file (parquet)",
	)]
        trajectories_path: String,

        #[arg(
	    short = 'm',
	    long = "mRNA_start_times",
	    help = "mRNA start times",
	    required = true
	)]
        mRNA_start_times_path: String,

        #[arg(
	    short = 'n',
	    long = "num_codons",
	    help = "Number of codons in the mRNA (must match value defined in \
		    simulation!)",
	    required = true
	)]
        num_codons: usize,

        #[arg(
	    short = 'o',
	    long = "out_dir",
	    default_value = "./",
	    help = "Output directory",
	    required = false
	)]
        out_dir: Option<String>,
	
	#[arg(short, long, action = ArgAction::Count)]
	verbosity: u8,

	#[arg(
	    short = 'L',
	    long = "footprint-length",
	    default_value_t = 10,
	    help = "Size of a ribosome footprint, in codons"
	)]
        footprint_length: usize,

	#[arg(
	    short = 'A',
	    long = "a-pos-length",
	    default_value_t = 4,
	    help = "Location of A-site position, in codons"
	)]
        a_pos_length: usize,

	#[arg(
	    short = 'W',
	    long = "warmup-time",
	    default_value_t = 0.0,
	    help = "Only collect statistics after a \"warmup\" time"
	)]
	warmup_time: f32
    },

    #[command(
	about = "Simulate a pulse-chase simulation using a pre-existing \
		 translation simulation dataset."
    )]
    PulseChase {
        #[arg(
	    short = 't',
	    long = "trajectories_path",
	    help = "Input trajectories file (parquet)",
	    required = true
	)]
        trajectories_path: String,

        #[arg(
	    short = 'o',
	    long = "out_file",
	    help = "Output file (parquet)",
	    required = true
	)]
        out_file: String,

	// TODO: allow loading sequence from a FASTA file?
        #[arg(
	    short = 's',
	    long = "sequence",
	    default_value = "MTMITDSLAVVLQRRDWENPGVTQLNRLAAHPPFASWRNSEEARTDRPSQQ\
			     LRSLNGEWRFAWFPAPEAVPESWLECDLPEADTVVVPSNWQMHGYDAPIYT\
			     NVTYPITVNPPFVPTENPTGCYSLTFNVDESWLQEGQTRIIFDGVNSAFHL\
			     WCNGRWVGYGQDSRLPSEFDLSAFLRAGENRLAVMVLRWSDGSYLEDQDMW\
			     RMSGIFRDVSLLHKPTTQISDFHVATRFNDDFSRAVLEAEVQMCGELRDYL\
			     RVTVSLWQGETQVASGTAPFGGEIIDERGGYADRVTLRLNVENPKLWSAEI\
			     PNLYRAVVELHTADGTLIEAEACDVGFREVRIENGLLLLNGKPLLIRGVNR\
			     HEHHPLHGQVMDEQTMVQDILLMKQNNFNAVRCSHYPNHPLWYTLCDRYGL\
			     YVVDEANIETHGMVPMNRLTDDPRWLPAMSERVTRMVQRDRNHPSVIIWSL\
			     GNESGHGANHDALYRWIKSVDPSRPVQYEGGGADTTATDIICPMYARVDED\
			     QPFPAVPKWSIKKWLSLPGETRPLILCEYAHAMGNSLGGFAKYWQAFRQYP\
			     RLQGGFVWDWVDQSLIKYDENGNPWSAYGGDFGDTPNDRQFCMNGLVFADR\
			     TPHPALTEAKHQQQFFQFRLSGQTIEVTSEYLFRHSDNELLHWMVALDGKP\
			     LASGEVPLDVAPQGKQLIELPELPQPESAGQLWLTVRVVQPNATAWSEAGH\
			     ISAWQQWRLAENLSVTLPAASHAIPHLTTSEMDFCIELGNKRWQFNRQSGF\
			     LSQMWIGDKKQLLTPLRDQFTRAPLDNDIGVSEATRIDPNAWVERWKAAGH\
			     YQAEAALLQCTADTLADAVLITTAHAWQHQGKTLFISRKTYRIDGSGQMAI\
			     TVDVEVASDTPHPARIGLNCQLAQVAERVNWLGLGPQENYPDRLTAACFDR\
			     WDLPLSDMYTPYVFPSENGLRCGTRELNYGPHQWRGDFQFNISRYSQQQLM\
			     ETSHRHLLHAEEGTWLNIDGFHMGIGGDDSWSPSVSAEFQLSAGRYHYQLV\
			     WCQK*",
	    help = "Amino acid sequence. Defaults to E. coli lacZ. Must include \
		    stop codon (*).",
	    conflicts_with = "fasta_path",
	)]
        sequence: Option<String>,

	#[arg(
	    short = 'f',
	    long = "fasta",
	    help = "Fasta file path, with amino acid sequence. Only use first sequence.",
	    conflicts_with = "sequence"
	)]
	fasta_path: Option<String>,

        #[arg(
	    short = 'l',
	    long = "labels",
	    default_value = "M",
	    help = "Radioactively labeled amino acid. Defaults to methionine.",
	)]
        labels: String,

        #[arg(
	    short = 'p',
	    long = "pulse_time",
	    help = "Pulse time. It makes sense to allow the polysome to populate \
		    before pulse-chase.",
	    required = true
	)]
        pulse_time: f32,

        #[arg(
	    short = 'c',
	    long = "chase_time",
	    help = "Chase time. Pulse duration should be short (e.g. 10 sec).",
	    required = true
	)]
        chase_time: f32,
    },
}

/*
   Boilerplate function that starts the program and immediately calls run(),
   and takes care of error handling.
*/
fn main() -> Result<()> {
    if let Err(e) = run() {
        eprintln!("Error: {:?}", e);
        std::process::exit(1);
    }
    Ok(())
}

/*
   Process command-line arguments & invoke the corresponding subroutine.
*/
fn run() -> Result<()> {
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Simulate {
	    simulation_time, out_dir, mRNA_half_life, num_mRNA, mRNA_freq,
	    rates_string, rates_file, events_string, events_file, num_threads,
	    verbosity, footprint_length, a_pos_length
	} => {
	    run_simulation(
		simulation_time,
		out_dir,
		mRNA_half_life,
		num_mRNA,
		mRNA_freq,
		rates_string,
		rates_file,
		events_string,
		events_file,
		num_threads,
		verbosity,
		footprint_length,
		a_pos_length
	    )
	}

        Commands::Statistics {
	    trajectories_path, mRNA_start_times_path, num_codons, out_dir,
	    verbosity, footprint_length, a_pos_length, warmup_time
	} => {
	    run_statistics(
		trajectories_path,
		mRNA_start_times_path,
		num_codons,
		out_dir,
		verbosity,
		footprint_length,
		a_pos_length,
		warmup_time
	    )
	}

        Commands::PulseChase {
	    trajectories_path, out_file, sequence, fasta_path, labels,
	    pulse_time, chase_time
	} => {
	    run_pulse_chase(
		trajectories_path,
		out_file,
		sequence,
		fasta_path,
		labels,
		pulse_time,
		chase_time,
	    )
	}
    }
}

fn run_pulse_chase(
    trajectories_path : String,
    out_file : String,
    sequence : Option<String>,
    fasta_path: Option<String>,
    labels : String,
    pulse_time : f32,
    chase_time : f32,
) -> Result<()> {

    let mut seq: String;
    match fasta_path {
	Some(fasta_path) => {
	    seq = util::read_first_protein_sequence(&fasta_path)?;
	}
	None => {
	    seq = sequence.unwrap();
	}
    }

    // Standardize on uppercase to simplify things
    seq = seq.to_uppercase();
    
    // Make sure that sequence ends with a stop codon, and only has 1 stop codon.
    if ! (seq.ends_with("*") && seq.matches("*").count() == 1) {
	bail!("Invalid stop codons in sequence");
    }
    
    // Read in the trajectories
    let r = File::open(&trajectories_path)?;
    let trajectories_df = ParquetReader::new(r).finish()?;

    let labels_vec : Vec<char> = labels
	.to_uppercase()
	.chars()
	.collect();

    let num_codons = seq.len() as u16;

    let out_path = PathBuf::from(&out_file);
    let of = File::create(out_path)?;

    // Calculate label positions based on sequence & labels
    let label_positions: Vec<u16> = seq
	.chars()
	.enumerate()
	.filter(|&(_, ch)| labels_vec.contains(&ch))
	.map(|(ind, _)| ind as u16)
	.collect();

    let mut pulse_chase = dataframe_operations::calculate_pulse_chase(
	trajectories_df, pulse_time, chase_time, label_positions, num_codons
    )?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut pulse_chase)?;

    Ok(())
}

// Does a directory exist? Is it writeable?
fn ensure_writable_dir<P: AsRef<Path>>(path: P, temp_filename : String) -> Result<()> {
    let path = path.as_ref();

    if path.exists() {
        if !path.is_dir() {
            anyhow::bail!("Path exists but is not a directory: {}", path.display());
        }

        // Check writability by trying to create a temp file
        let test_path = path.join(temp_filename);
        write(&test_path, b"test").context("Directory exists but is not writable")?;
        remove_file(&test_path).ok(); // Cleanup, ignore error
    } else {
        create_dir_all(&path)
            .with_context(|| format!("Failed to create directory: {}", path.display()))?;
    }

    Ok(())
}


fn run_statistics(
    trajectories_path : String,
    mRNA_start_times_path : String,
    num_codons : usize,
    out_dir : Option<String>,
    verbosity : u8,
    footprint_length: usize,
    a_pos_length: usize,
    warmup_time: f32
) -> Result<()> {
    // Did the user specify out_dir? If not, then default to current working directory:
    let out_dir = match out_dir {
	Some(out_dir) => out_dir,
	None => env::current_dir()
	    .expect("Failed to get current directory")
	    .to_string_lossy()
	    .into()
    };
    
    // Make sure that output directory exists and can be written to
    ensure_writable_dir(PathBuf::from(&out_dir), "step_times.pq".to_string())?;
    
    let num_codons = num_codons as u16;

    // Read in the trajectories
    let r = File::open(&trajectories_path)?;
    let trajectories_df = ParquetReader::new(r).finish()?;

    // Read in the mRNA start times
    let r = File::open(&mRNA_start_times_path)?;
    let mRNA_start_times_df = ParquetReader::new(r).finish()?;
    let num_mRNA = mRNA_start_times_df.height() as u16;

    // Step times
    if verbosity >= 1 {
	println!("Step times");
    }
    let mut step_times = dataframe_operations::calculate_complete_step_times(
	&trajectories_df, &mRNA_start_times_df, warmup_time
    )?;

    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("step_times.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut step_times)?;

    // Step times (mean)
    if verbosity >= 1 {
	println!("Mean step times");
    }
    // Read in complete step times
    let mut in_path = PathBuf::from(&out_dir);
    in_path.push("step_times.pq");
    let file_st = File::open(&in_path)?;
    let df_st = ParquetReader::new(file_st).finish()?;
    let mut mean_step_times = dataframe_operations::calculate_mean_step_times(
	&df_st
    )?;
    
    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("mean_step_times.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut mean_step_times)?;

    // Occupancy
    if verbosity >= 1 {
	println!("Final occupancy");
    }
    let mut occupancy = dataframe_operations::calculate_final_occupancy(
	&trajectories_df, num_codons, num_mRNA, footprint_length, a_pos_length
    )?;

    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("occupancy.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut occupancy)?;


    // Protein Production
    
    // 1. By mRNA
    if verbosity >= 1 {
	println!("Protein production by mRNA");
    }
    let mut prot_mRNA = dataframe_operations::calculate_protein_production_by_mrna(
	&trajectories_df, num_codons, warmup_time
    )?;

    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("protein_production_by_mRNA.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut prot_mRNA)?;

    // 2. Total
    if verbosity >= 1 {
	println!("Protein production");
    }
    let mut total_mRNA = dataframe_operations::calculate_total_protein_production(
	&trajectories_df, num_codons, warmup_time
    )?;

    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("protein_production.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut total_mRNA)?;

    // Synthesis times
    if verbosity >= 1 {
	println!("Synthesis times");
    }
    let mut synthesis = dataframe_operations::calculate_synthesis_times(
	&trajectories_df, num_codons, warmup_time
    )?;

    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("synthesis_times.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut synthesis)?;

    // Collisions
    if verbosity >= 1 {
	println!("Collisions");
    }
    let mut collisions = dataframe_operations::calculate_collision_frequency(
	&trajectories_df, num_codons, footprint_length, warmup_time
    )?;

    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("collisions.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut collisions)?;

    // Nascent peptide length distribution
    if verbosity >= 1 {
	println!("Nascent peptide length distribution");
    }
    let mut lengths = dataframe_operations::calculate_final_nascent_peptide_length_distribution(
	&trajectories_df, num_codons
    )?;

    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("length_distribution.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut lengths)?;

    // Ribosome spacing
    if verbosity >= 1 {
	println!("Ribosome spacing");
    }
    let mut spacings = dataframe_operations::calculate_final_average_ribosome_spacing(
	&trajectories_df, num_codons
    )?;

    let mut out_path = PathBuf::from(&out_dir);
    out_path.push("spacing.pq");
    let of = File::create(out_path)?;

    let writer = ParquetWriter::new(of);
    writer.finish(&mut spacings)?;
    
    Ok(())
}

fn run_simulation(
    simulation_time : f32,
    out_dir : Option<String>,
    mRNA_half_life : f32,
    num_mRNA : Option<u16>,
    mRNA_freq : Option<f32>,
    rates_string : Option<String>,
    rates_file : Option<String>,
    events_string : Option<String>,
    events_file : Option<String>,
    num_threads : usize,
    verbosity : u8,
    footprint_length: usize,
    a_pos_length: usize
) -> Result<()> {
    // Sanity checks
    assert!(footprint_length > a_pos_length);
    // A, P, E. Nonsensical to have a footprint less than 3 codons!
    assert!(a_pos_length >= 3);
    
    // Parse translation rates either from a string or from a file
    let rates: Vec<f32>;
    if let Some(rates_string) = rates_string {
	rates = util::parse_rates(
	    &rates_string
	)?;
    } else if let Some(rates_file) = rates_file {
	rates = util::init_rates_from_file(
	    &rates_file
	)?;
    } else {
	bail!("no probabilities specified");
    }

    /* Make sure there are at least 3 rates: 1 start codon, 1 regular codon, 1 
       stop codon. Otherwise it's not possible to undergo a complete translation
       cycle. */
    if rates.len() < 3 {
	bail!("Not enough codons!");
    }

    // Scheduled events:
    // At a given time, change the rate at a given position to a new value
    let scheduled_events : Vec<_>;
    if let Some(events_file) = events_file {
	scheduled_events = util::init_events_tsv(&events_file)
	    .context("Failed to read events file")?;
    } else if let Some(events_string) = events_string {
	scheduled_events = util::parse_events_str(&events_string)
	    .context("Failed to parse events string")?;
    }
    else {
	scheduled_events = Vec::new();
    }
    
    // Convert mRNA half-life to decay rate constant
    let mRNA_decay_rate: f32;
    if mRNA_half_life == f32::INFINITY {
	mRNA_decay_rate = 0.0_f32;
    } else {
	mRNA_decay_rate = 2.0_f32.ln() / mRNA_half_life;
    }

    // Did the user specify out_dir? If not, then default to current working
    // directory:
    let out_path : PathBuf;
    if let Some(out_dir) = out_dir {
	out_path = PathBuf::from(&out_dir);
    } else {
	out_path = env::current_dir()
	    .expect("Failed to get current directory");
    }
    
    // Make sure that output directory exists and can be written to
    ensure_writable_dir(&out_path, "trajectories.pq".to_string())?;
    
    /* This channel is used to transfer information in & out of the
       threadpool. */
    let (tx_mRNA, rx_mRNA) = channel();

    let pool = ThreadPool::new(num_threads);

    // Use this to record when each mRNA starts
    let mut mRNA_start_times_record: Vec<f32> = Vec::new();

    // Induce
    if let Some(mRNA_freq) = mRNA_freq {

	if verbosity >= 1 {
	    println!("Spawning mRNA with frequency {freq:.*}", 2, freq=mRNA_freq);
	}
	
	let mut current_mRNA = 0;
	let mut rng = rand::rng();
	let mRNA_spawn_pdf = Exp::new(mRNA_freq)?;

	let mut mRNA_time = mRNA_spawn_pdf.sample(&mut rng);

	while mRNA_time < simulation_time {
	    // Each mRNA starts at a different time
	    mRNA_start_times_record.push(mRNA_time);
	    
	    let tx_mRNA = tx_mRNA.clone();
	    let mut rates = rates.clone();
	    let mut scheduled_events = scheduled_events.clone();

	    /* call run simplton here */
	    pool.execute(move || {
		let trajectories = simplton::run_simplton(
		    mRNA_time,
		    simulation_time,
		    &mut rates,
		    current_mRNA,
		    mRNA_decay_rate,
		    &mut scheduled_events,
		    footprint_length,
		    a_pos_length
		);

		let return_message = SimulationMessage {
		    mRNA_number: current_mRNA,
		    trajectories: trajectories,
		};

		tx_mRNA
		    .send(return_message)
		    .expect("Message channel should be open");
	    });

	    // Advance mRNA time to the next spawned mRNA
	    mRNA_time = mRNA_spawn_pdf.sample(&mut rng) + mRNA_time;
	    current_mRNA += 1;
	}
    }
    // Constant # mRNA
    else if let Some(num_mRNA) = num_mRNA {
	if verbosity >= 1 {
	    println!("Generating {} mRNAs", num_mRNA);
	}
	
	for current_mRNA in 0..num_mRNA {
	    // All mRNAs start at time zero
	    mRNA_start_times_record.push(0.0);

	    let tx_mRNA = tx_mRNA.clone();
	    let mut rates = rates.clone();
	    let mut scheduled_events = scheduled_events.clone();

	    pool.execute(move || {
		/* call run simplton here */
		let trajectories = simplton::run_simplton(
		    0.0,
		    simulation_time,
		    &mut rates,
		    current_mRNA,
		    mRNA_decay_rate,
		    &mut scheduled_events,
		    footprint_length,
		    a_pos_length
		);

		let return_message = SimulationMessage {
		    mRNA_number: current_mRNA,
		    trajectories: trajectories,
		};

		tx_mRNA
		    .send(return_message)
		    .expect("Message channel should be open");
	    });
	    
	}
    }

    drop(tx_mRNA);

    /* Parquet file creation for ribosome trajectories */
    let num_columns = 4;
    let mut schema_vec: Vec<Field> = Vec::with_capacity(num_columns);
    schema_vec.push(Field::new("mRNA", DataType::UInt32, false));
    schema_vec.push(Field::new("ribosome", DataType::UInt32, false));
    schema_vec.push(Field::new("pos", DataType::UInt16, false));
    schema_vec.push(Field::new("timestamp", DataType::Float32, false));

    // Create the parquet file
    let schema = Schema::from(schema_vec);

    let options = WriteOptions {
	write_statistics: true,
	compression: CompressionOptions::Snappy,
	version: Version::V2,
	data_pagesize_limit: None,
    };

    let out_file = BufWriter::new(File::create(
	PathBuf::from(&out_path).join("trajectories.pq"),
    )?);

    let mut writer = FileWriter::try_new(out_file, schema.clone(), options)?;

    /* Done with trajectories parquet file creation */

    /* The encodings */
    /* https://docs.rs/arrow2/0.18.0/arrow2/io/parquet/write/enum.Encoding.html
       NOTE: Since order of results isn't guaranteed, mRNA # & ribosome & mightn't
       be in sorted order. --> Use Plain for first 3 columns.
       Or, would DeltaBinaryPacked work well? */
    let encodings_vec = vec![
	vec![Encoding::Plain], // mRNA #
	vec![Encoding::Plain], // ribosome #
	vec![Encoding::Plain], // pos
	vec![Encoding::Plain], // timestamp
    ];

    /* Process each message: a completed simulation */
    let mut terminal = term::stdout().context("Failed to create stdout terminal")?;
    for message in rx_mRNA.iter() {
	terminal.carriage_return().context("Failed carriage return")?;
	terminal.delete_line().context("Failed delete line")?;

	// Handle errors on trajectories here, in the main thread
	let trajectories = message.trajectories?;
	if verbosity >= 1 {
	    write!(terminal, "mRNA # {}", message.mRNA_number)?;
	}

	// Write the incoming trajectories to the parquet file
	let row_groups = RowGroupIterator::try_new(
	    vec![Ok(trajectories)].into_iter(),
	    &schema,
	    options,
	    encodings_vec.clone(),
	)?;

	for group in row_groups {
	    writer.write(group?)?;
	}

	io::stdout().flush().context("Failed to flush stdout")?;
    }

    // Close the trajectories parquet file
    let _ = writer.end(None)?;

    // Write the mRNA start times to parquet
    write_mRNA_start_times(&out_path, mRNA_start_times_record)?;

    /* Simulations complete */
    if verbosity >= 1 {
	writeln!(terminal, "")?;
    }

    Ok(())
}


/* Parquet file creation for mRNA start times */
pub fn write_mRNA_start_times(
    out_path: &PathBuf,
    mRNA_start_times_record: Vec<f32>
) -> Result<()> {
    /*
    Parquet file creation for mRNA start times
    Needed for some statistics, e.g. initiation time of pioneer ribosome,
    pulse-chase simulation, protein production.
    - If fixed number of mRNA, then will always be 0.0.
    - If spawned, then will be variable.

    Column 1: mRNA #
    Column 2: start time
    */

    let options = WriteOptions {
        write_statistics: true,
        compression: CompressionOptions::Snappy,
        version: Version::V2,
        data_pagesize_limit: None,
    };

    let mRNA_num_columns = 2;
    let mut mRNA_schema_vec: Vec<Field> = Vec::with_capacity(mRNA_num_columns);

    mRNA_schema_vec.push(Field::new("mRNA", DataType::UInt32, false));
    mRNA_schema_vec.push(Field::new("start time", DataType::Float32, false));

    let schema_mRNA = Schema::from(mRNA_schema_vec);

    let out_file_mRNA = BufWriter::new(File::create(
        PathBuf::from(&out_path).join("mRNA_start_times.pq"),
    )?);

    let mut mRNA_writer = FileWriter::try_new(
	out_file_mRNA, schema_mRNA.clone(), options
    )?;

    // Return a "chunk" of arrays
    let mut arrays_mRNA: Vec<Arc<dyn Array>> = Vec::with_capacity(
	mRNA_num_columns
    );

    // Give every ribosome entry a copy of its mRNA simulation number, too.
    let num_mRNA = mRNA_start_times_record.len();
    let mRNA_number_list: Vec<u16> = (0..num_mRNA as u16).collect();
    let mRNA_number_array = PrimitiveArray::from_vec(mRNA_number_list);
    arrays_mRNA.push(mRNA_number_array.arced());

    let mRNA_start_times_array = PrimitiveArray::from_vec(
	mRNA_start_times_record
    );

    arrays_mRNA.push(mRNA_start_times_array.arced());

    let mRNAs_chunk = Chunk::try_new(arrays_mRNA);

    let encodings_vec_mRNA = vec![
        vec![Encoding::Plain], // mRNA #
        vec![Encoding::Plain], // timestamp
    ];

    // Write the incoming trajectories to the parquet file
    let row_groups = RowGroupIterator::try_new(
        vec![mRNAs_chunk].into_iter(),
        &schema_mRNA,
        options,
        encodings_vec_mRNA.clone(),
    )?;

    for group in row_groups {
        mRNA_writer.write(group?)?;
    }

    // Close the parquet file
    let _ = mRNA_writer.end(None)?;

    Ok(())
}
