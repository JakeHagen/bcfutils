mod commands;
use commands::*;
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
        #[clap(long, value_parser, default_value_t = 1)]
        threads: usize,
    },
    MNV {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
    },
    GLKUP {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
        #[clap(long, short)]
        fields: Option<String>,
        #[clap(long, takes_value = false)]
        dbnsfp: bool,
        #[clap(long, value_parser, default_value_t = 1)]
        threads: usize,
    },
    QC {
        input: Option<String>,
        #[clap(long, short)]
        output: Option<String>,
        #[clap(long, short)]
        field: Option<String>,
        #[clap(long, short)]
        rare: f32,
        #[clap(long, short)]
        gq: i32,
        //#[clap(long, short)]
        //bed: Option<String>,
        #[clap(long, value_parser, default_value_t = 1)]
        threads: usize,
    }
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
        Commands::MCSQ { input, output, gff, threads } => {
            mcsq::mcsq(input.as_deref(), output.as_deref(), gff.as_deref(), threads)
        }
        Commands::MNV { input, output } => {
            mnv::mnv(input.as_deref(), output.as_deref())
        }
        Commands::GLKUP { input, output, fields, dbnsfp, threads } => {
            glkup::glkup(input.as_deref(), output.as_deref(), fields.as_deref(), dbnsfp, threads)
        },
        Commands::QC { input, output, field, rare, gq, threads } => {
            qc::qc(input.as_deref(), output.as_deref(), field.as_deref(), rare, gq, threads)
        }
    }
}
