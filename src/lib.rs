extern crate rust_htslib;
use crate::rust_htslib::bcf::{Reader, Format, Writer};

pub fn get_rdr(input: Option<&str>) -> rust_htslib::bcf::Reader {
    match input {
        None => {
           return  Reader::from_stdin().expect("failed to create reader from stdin");
        }
        Some("-") => {
            return Reader::from_stdin().expect("failed to create reader from stdin");
        }
        Some(inner) => {
            return Reader::from_path(inner).expect("failed to create reader from file");
        }
    }
}

pub fn get_wrtr(input: Option<&str>, hdr: &rust_htslib::bcf::Header) -> rust_htslib::bcf::Writer {
    match input {
        None => {
            return Writer::from_stdout(hdr, false, Format::Bcf).expect("unable to create writer on stdout");
        }
        Some("-") => {
            return Writer::from_stdout(hdr, false, Format::Bcf).expect("unable to create writer on stdout");
        }
        Some(inner) => {
            return Writer::from_path(inner, hdr, false, Format::Bcf).expect("unable to create writer with filepath");
        }
    }
}


