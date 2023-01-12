use linear_map::LinearMap;
use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::record::Buffer;
use rust_htslib::bcf::{Header, Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::str;

fn get_field<'a>(s: &'a str, field: &'a str) -> Option<&'a str> {
    let sb = match s.find(field) {
        Some(b) => b + field.len() + 1, // +1 for the equal sign
        None => return None,
    };
    let eb = match s[sb..].find(";") {
        Some(b) => sb + b,
        None => s.len(),
    };
    return Some(&s[sb..eb]);
}

fn has_tag(tags: &str, tag: &str) -> bool {
    match tags.find(tag) {
        Some(_) => return true,
        None => return false,
    };
}

fn extract_appris(tags: &str) -> Option<&str> {
    let sb = match tags.find("appris_") {
        Some(b) => b + 7, // 7 == length of appris_
        None => return None,
    };
    let eb = match tags[sb..].find(",") {
        Some(b) => sb + b,
        None => tags.len(),
    };
    return Some(&tags[sb..eb]);
}

fn extract_uncertain_start_end(tags: &str) -> String {
    let mut found_tags = vec![];
    let possible_tags = vec!["cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"];
    for tag in possible_tags {
        if has_tag(tags, tag) {
            found_tags.push(tag);
        }
    }
    return found_tags.join("&");
}

fn add_trns_to_map(line: String, map: &mut HashMap<String, String>) {
    for (i, col) in line.split("\t").enumerate() {
        if i == 2 {
            if col != "transcript" {
                return;
            }
        }
        if i == 8 {
            let mut s = String::new();
            match get_field(col, "gene_id") {
                Some(gene_id) => {
                    s.push('|');
                    s.push_str(gene_id);
                }
                None => s.push('|'),
            };
            match get_field(col, "tag") {
                Some(tags) => {
                    // add CANONICAL
                    if has_tag(tags, "Ensembl_canonical") {
                        s.push_str("|YES")
                    } else {
                        s.push('|')
                    }
                    // add appris
                    match extract_appris(tags) {
                        Some(appris) => {
                            s.push('|');
                            s.push_str(appris);
                        }
                        None => s.push('|'),
                    }
                    // add CCDS
                    if has_tag(tags, "CCDS") {
                        s.push_str("|CCDS")
                    } else {
                        s.push('|')
                    }
                    //add readthrough
                    if has_tag(tags, "readthrough_transcript") {
                        s.push_str("|readthrough_transcript")
                    } else {
                        s.push('|')
                    }
                    // add uncertain start/end
                    s.push('|');
                    s.push_str(&extract_uncertain_start_end(tags));
                }
                None => s.push_str("|||||"),
            };
            // add TSL
            match get_field(col, "transcript_support_level") {
                Some(tsl) => {
                    s.push('|');
                    s.push_str(tsl);
                }
                None => s.push('|'),
            };

            match get_field(col, "transcript_id") {
                Some(trns_id) => {
                    s.push('|');
                    s.push_str(trns_id);
                    let eb = trns_id
                        .find(".")
                        .expect("transcript doesnt have '.' notation");
                    map.insert(trns_id[0..eb].to_string(), s);
                }
                None => s.push('|'),
            }
        }
    }
}

fn build_trx_map(gff_fp: Option<&str>) -> HashMap<String, String> {
    let gff = File::open(gff_fp.unwrap());
    let rdr = BufReader::new(gff.expect("couldnt make buffer for gff file"));

    let mut trns_map: HashMap<String, String> = HashMap::new();

    for line in rdr.lines() {
        if let Ok(l) = line {
            if &l[0..1] == "#" {
                continue;
            }
            add_trns_to_map(l, &mut trns_map);
        }
    }
    return trns_map;
}

fn get_bcsq_hdr_map(hdr_recs: Vec<HeaderRecord>) -> Option<LinearMap<String, String>> {
    for hrec in hdr_recs.iter() {
        match hrec {
            HeaderRecord::Info { key, values } => {
                if values.get("ID").unwrap() == "BCSQ" {
                    return Some(values.clone());
                }
            }
            _ => continue,
        }
    }
    return None;
}

fn get_num_bcsq_keys(desc: &str) -> usize {
    return desc.split('|').count();
}

fn get_canon_rank(csq: &Vec<&str>) -> u8 {
    let canon = csq[8];
    return match canon {
        "YES" => 1,
        _ => 99,
    };
}

fn get_appris_rank(csq: &Vec<&str>) -> u8 {
    let appris = csq[9];
    return match appris {
        "principal_1" => 1,
        "principal_2" => 2,
        "principal_3" => 3,
        "principal_4" => 4,
        "principal_5" => 5,
        "alternative_1" => 6,
        "alternative_2" => 7,
        _ => 99,
    };
}

fn get_tsl_rank(csq: &Vec<&str>) -> u8 {
    let tsl = csq[13];
    return match tsl {
        "1" => 1,
        "2" => 2,
        "3" => 3,
        "4" => 4,
        "5" => 5,
        _ => 99,
    };
}

fn get_severity_rank(csq: &Vec<&str>) -> u8 {
    let mut max: u8 = 99;
    for c in csq[0].split("&") {
        let rank = match c {
            "transcript_ablation" => 1,
            "splice_acceptor" => 2,
            "splice_donor" => 2,
            "stop_gained" => 3,
            "frameshift" => 3,
            "stop_lost" => 3,
            "start_lost" => 3,
            "disruptive" => 4,
            "exon_loss" => 5,
            "transcript_amplification" => 6,
            "inframe_altering" => 7,
            "inframe_insertion" => 7,
            "inframe_deletion" => 7,
            "missense" => 7,
            "protein_altering" => 7,
            "inframe" => 7,
            "splice_region" => 8,
            "incomplete_terminal_codon" => 9,
            "synonymous" => 10,
            "stop_retained" => 10,
            "start_retained" => 10,
            "coding_sequence" => 11,
            "mature_miRNA" => 11,
            "5_prime_utr" => 12,
            "3_prime_utr" => 12,
            "non_coding_transcript_exon" => 13,
            "intron" => 14,
            "NMD_transcript" => 14,
            "non_coding" => 15,
            "non_coding_transcript" => 15,
            "downstream" => 16,
            "upstream" => 16,
            "TF_binding_site" => 17,
            "TFBS" => 17,
            "regulatory" => 18,
            "feature_truncation" => 19,
            "feature_elongation" => 19,
            "intergenic" => 20,
            _ => 99,
        };
        if rank < max {
            max = rank
        }
    }
    return max;
}

fn get_biotype_rank(csq: &Vec<&str>) -> u8 {
    return match csq[3] {
        "protein_coding" => 1,
        _ => 99,
    };
}

fn get_readthrough_rank(csq: &Vec<&str>) -> u8 {
    return match csq[11] {
        "readthrough_transcript" => 99,
        _ => 1,
    };
}

const CANON: usize = 0;
const APPRIS: usize = 1;
const TSL: usize = 2;
const SEVERE: usize = 3;
const BIOTYPE: usize = 4;
const READTHROUGH: usize = 5;

fn get_ranks(csq: &Vec<&str>) -> [u8; 6] {
    return [
        get_canon_rank(csq),
        get_appris_rank(csq),
        get_tsl_rank(csq),
        get_severity_rank(csq),
        get_biotype_rank(csq),
        get_readthrough_rank(csq),
    ];
}

fn compare_ranks(
    p_rank: &mut [u8; 6],
    p_csq_idx: &mut usize,
    c_rank: [u8; 6],
    c_csq_idx: usize,
    comps: Vec<usize>,
) {
    for comp in comps {


        if c_rank[comp] < p_rank[comp] {
            *p_rank = c_rank;
            *p_csq_idx = c_csq_idx;
            return;
        }
        if c_rank[comp] > p_rank[comp] {
            return;
        }
    }
}

pub fn mcsq(input: Option<&str>, output: Option<&str>, gff_fp: Option<&str>, threads: &usize) {
    let mut bcf = bcfutils::get_rdr(input);
    bcf.set_threads(threads.clone())
        .expect("unable to set reader threads");

    let hdrv = bcf.header();

    let bcsq_map =
        get_bcsq_hdr_map(hdrv.header_records()).expect("was not able to get BCSQ header info line");
    let num_keys = get_num_bcsq_keys(bcsq_map.get("Description").expect("bcsq map doesnt have \"Description\", it really should though, something funky is happening"));

    let mut hdr = Header::from_template(&hdrv);

    hdr.remove_info(b"BCSQ");
    hdr.push_record(r#"##INFO=<ID=BCSQ,Number=.,Type=String,Description="Local consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.  html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change|gene_id|CANONICAL|appris|ccds|unknown_start_end|TSL|transcript_id">"#.as_bytes());

    let bcsq_fields = vec![
        "Consequence",
        "gene",
        "transcript",
        "biotype",
        "strand",
        "amino_acid_change",
        "dna_change",
        "gene_id",
        "CANONICAL",
        "appris",
        "ccds",
        "readthrough",
        "unknown_start_end",
        "TSL",
        "transcript_id",
    ];

    for new_field in &bcsq_fields {
        hdr.push_record(format!("##INFO=<ID=canon_{},Number=1,Type={},Description=\"canon {}\">", new_field, "String", new_field).as_bytes());
        hdr.push_record(format!("##INFO=<ID=pick_{},Number=1,Type={},Description=\"picked csq, canonical->appris->TSL->biotype->severity {}\">", new_field, "String", new_field).as_bytes());
        hdr.push_record(format!("##INFO=<ID=worst_{},Number=1,Type={},Description=\"worst {}\">", new_field, "String", new_field).as_bytes());
        hdr.push_record(format!("##INFO=<ID=wpc_{},Number=1,Type={},Description=\"worst protein coding {}\">", new_field, "String", new_field).as_bytes());
    }

    let mut obcf = bcfutils::get_wrtr(output, &hdr);
    obcf.set_threads(threads.clone())
        .expect("unable to set writer threads");

    let trx_map = build_trx_map(gff_fp);
    let mut b = Buffer::new();
    for record_result in bcf.records() {
        let mut record = record_result.expect("fail to read record");
        obcf.translate(&mut record);
        let bcsqs = match record.info_shared_buffer(b"BCSQ", &mut b).string().unwrap() {
            Some(b) => b,
            None => {
                obcf.write(&record).expect("failed to write record");
                continue;
            }
        };

        let mut mcsqs = vec![];

        let mut p_csq_idx = 0;
        let mut p_rank = [99, 99, 99, 99, 99, 99];
        let mut w_csq_idx = 0;
        let mut w_rank = [99, 99, 99, 99, 99, 99];
        let mut wpc_csq_idx = 0;
        let mut wpc_rank = [99, 99, 99, 99, 99, 99];

        for (i, bcsq_b) in bcsqs.iter().enumerate() {
            let mut mcsq = String::new();
            let bcsq = str::from_utf8(bcsq_b).unwrap();
            let num_c_bcsq = get_num_bcsq_keys(bcsq);

            mcsq.push_str(bcsq);

            let num_to_fill = num_keys - num_c_bcsq;
            for _ in 0..num_to_fill {
                mcsq.push('|');
            }
            let trn = bcsq.split('|').collect::<Vec<&str>>()[2];
            match trx_map.get(trn) {
                Some(t) => {
                    mcsq.push_str(t);
                }
                None => mcsq.push_str("||||||||"),
            }

            let r = get_ranks(&mcsq.split("|").collect::<Vec<&str>>());
            mcsqs.push(mcsq);

            compare_ranks(
                &mut p_rank,
                &mut p_csq_idx,
                r,
                i,
                vec![READTHROUGH, CANON, APPRIS, TSL, BIOTYPE, SEVERE],
            );

            compare_ranks(
                &mut w_rank,
                &mut w_csq_idx,
                r,
                i,
                vec![READTHROUGH, SEVERE, CANON, APPRIS, TSL, BIOTYPE],
            );

            compare_ranks(
                &mut wpc_rank,
                &mut wpc_csq_idx,
                r,
                i,
                vec![READTHROUGH, BIOTYPE, SEVERE, CANON, APPRIS, TSL],
            );
        }
        
        for (i, f) in mcsqs[p_csq_idx].split("|").enumerate() {
            record
                .push_info_string(
                    format!("pick_{}", bcsq_fields[i]).as_bytes(),
                    &[f.as_bytes()],
                )
                .expect("failed to set canon_BCSQ field");
        }
        if p_rank[0] == 1 {
            for (i, f) in mcsqs[p_csq_idx].split("|").enumerate() {
                record
                    .push_info_string(
                        format!("canon_{}", bcsq_fields[i]).as_bytes(),
                        &[f.as_bytes()],
                    )
                    .expect("failed to set canon_BCSQ field");
            }
        }
        for (i, f) in mcsqs[w_csq_idx].split("|").enumerate() {
            record
                .push_info_string(
                    format!("worst_{}", bcsq_fields[i]).as_bytes(),
                    &[f.as_bytes()],
                )
                .expect("failed to set canon_BCSQ field");
        }
        for (i, f) in mcsqs[wpc_csq_idx].split("|").enumerate() {
            record
                .push_info_string(
                    format!("wpc_{}", bcsq_fields[i]).as_bytes(),
                    &[f.as_bytes()],
                )
                .expect("failed to set canon_BCSQ field");
        }
      
        record
            .push_info_string(b"BCSQ", &[mcsqs.join(",").as_bytes()])
            .expect("failed to set QD info field");
        obcf.write(&record).expect("failed to write record");
    }
}
