extern crate rust_htslib;
extern crate linear_map;
mod commands;

use commands::*;
use std::str;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    QD {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
    },
    FamFreq {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
        #[clap(long, short)]
        pedigree: Option<String>,
    },
    MCSQ {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
        #[clap(long, short)]
        gff: Option<String>,
    },
    PICK {
          input: Option<String>,
          #[clap(long, short)]
          output: Option<String>,
      },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::QD { input, output } => {
            ann_qd::ann_qd(input.as_deref(), output.as_deref())
        }
        Commands::FamFreq { input, output, pedigree } => {
            fam_freq::fam_freq(pedigree.as_deref())
        }
        Commands::MCSQ { input, output, gff } => {
            mcsq::mcsq(input.as_deref(), output.as_deref(), gff.as_deref())
        }
        Commands::PICK { input, output } => {
            pick::pick(input.as_deref(), output.as_deref())
        }
    }
}
