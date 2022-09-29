use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::str;
use rust_htslib::bcf::{Read, Header};
use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::record::Buffer;
use serde::Deserialize;
use std::error::Error;
use std::process;
use std::io;

static DBNSFP: &[u8] = include_bytes!("../../lookups/dbNSFP4.2_gene");
static GNOMAD: &[u8] = include_bytes!("../../lookups/gnomad.v2.1.1.lof_metrics.by_gene.txt");

#[derive(Debug, Deserialize)]
struct Grow {
    gene: String,
    #[serde(rename = "pLI", deserialize_with = "csv::invalid_option")]
    gnomAD_pLI: Option<f32>,
    #[serde(deserialize_with = "csv::invalid_option")]
    oe_lof: Option<f32>,
    #[serde(deserialize_with = "csv::invalid_option")]
    oe_lof_upper: Option<f32>,
    #[serde(deserialize_with = "csv::invalid_option")]
    syn_z: Option<f32>,
    #[serde(deserialize_with = "csv::invalid_option")]
    mis_z: Option<f32>,
    #[serde(deserialize_with = "csv::invalid_option")]
    lof_z: Option<f32>,
    #[serde(deserialize_with = "csv::invalid_option")]
    exac_pLI: Option<f32>,
}

#[derive(Debug, Deserialize)]
struct DBrow {
    #[serde(rename = "Gene_name")]
    gene:String,
    #[serde(rename = "Gene_other_names")]
    gene_syn: String,
    #[serde(rename = "Gene_full_name")]
    gene_full: String,
    #[serde(rename = "Pathway(Uniprot)")]
    pathway_uniprot: String,
    #[serde(rename = "Pathway(BioCarta)")]
    pathway_biocarta: String,
    #[serde(rename = "Pathway(ConsensusDB)")]
    pathway_consensusPathDB: String,
    #[serde(rename = "Pathway(KEGG)_full")]
    pathway_kegg: String,
    #[serde(rename = "Function_description")]
    gene_function: String,
    #[serde(rename = "Disease_description")]
    gene_disease: String,
    #[serde(rename = "MIM_phenotype_id")]
    MIM_phenotype_id: String,
    #[serde(rename = "MIM_disease")]
    MIM_disease: String,
    #[serde(rename = "Orphanet_disorder_id")]
    orphanet_id: String,
    #[serde(rename = "Orphanet_disorder")]
    orphanet_disorder: String,
    #[serde(rename = "Orphanet_association_type")]
    orphanet_assoc_type: String,
    #[serde(rename = "Trait_association(GWAS)")]
    GWAS_trait: String,
    #[serde(rename = "HPO_id")]
    HPO_id: String,
    #[serde(rename = "HPO_name")]
    HPO_name: String,
    #[serde(rename = "GO_biological_process")]
    GO_bio_process: String,
    #[serde(rename = "GO_cellular_component")]
    GO_cellular_comp: String,
    #[serde(rename = "GO_molecular_function")]
    GO_molecular_func: String,
    #[serde(rename = "Tissue_specificity(Uniprot)")]
    UNIPROT_tissue_specificity: String,
    #[serde(rename = "Expression(egenetics)")]
    egenetics_expression: String,
    #[serde(rename = "Expression(GNF/Atlas)")]
    GNF_atlas_expression: String,
    #[serde(rename = "MGI_mouse_gene")]
    MGI_mouse_gene: String,
    #[serde(rename = "MGI_mouse_phenotype")]
    MGI_mouse_phenotype: String,
}


fn build_gnomad_map() -> HashMap<String, Grow> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(GNOMAD);
    let mut map: HashMap<String, Grow> = HashMap::new();
    for result in rdr.deserialize() {
          match result {
              Ok::<Grow, _>(grow) => {
                  map.insert(grow.gene.clone(), grow);
              },
              Err(e) => continue,
          }
      }
    return map;
}

fn build_dbnsfp_map() -> HashMap<String, DBrow> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(DBNSFP);
    let mut map: HashMap<String, DBrow> = HashMap::new();
    for result in rdr.deserialize() {
        match result {
            Ok::<DBrow, _>(dbrow) => {
                map.insert(dbrow.gene.clone(), dbrow);
            },
            Err(e) => continue,
        }
    }
    return map;
}

fn add_gnomad_fields(grow: &Grow, record: &mut rust_htslib::bcf::record::Record) {
    record.push_info_string("gnomAD_gene".as_bytes(), &[grow.gene.as_bytes()]).expect("failed to set gnomAD_gene field");
    match grow.gnomAD_pLI {
        Some(gnomad_pli) => {
            record.push_info_float("gnomAD_pLI".as_bytes(), &[gnomad_pli]).expect("failed to set gnomAD_pLI field");
        },
        None => (),
    };

    match grow.oe_lof {
        Some(oe_lof) => {
            record.push_info_float("oe_lof".as_bytes(), &[oe_lof]).expect("failed to set oe_lof field");
        },
        None => (),
    };
    match grow.oe_lof_upper {
        Some(oe_lof_upper) => {
            record.push_info_float("oe_lof_upper".as_bytes(), &[oe_lof_upper]).expect("failed to set oe_lof_upper field");
        },
        None => (),
    };
    match grow.syn_z {
        Some(syn_z) => {
            record.push_info_float("syn_z".as_bytes(), &[syn_z]).expect("failed to set syn_z field");
        },
        None => (),
    };
    match grow.mis_z {
        Some(mis_z) => {
            record.push_info_float("mis_z".as_bytes(), &[mis_z]).expect("failed to set mis_z field");
        },
        None => (),
    };
    match grow.lof_z {
        Some(lof_z) => {
            record.push_info_float("lof_z".as_bytes(), &[lof_z]).expect("failed to set lof_z field");
        },
        None => (),
    };
    match grow.exac_pLI {
        Some(exac_pli) => {
            record.push_info_float("exac_pLI".as_bytes(), &[exac_pli]).expect("failed to set exac_pLI field");
        },
        None => (),
    };
    return
}

fn add_dbnsfp_fields(dbrow: &DBrow, record: &mut rust_htslib::bcf::record::Record) {
    record.push_info_string("dbnsfp_gene".as_bytes(), &[dbrow.gene.as_bytes()]).expect("failed to set dbnsfp_gene field");
    if dbrow.gene_syn != "." {
        record.push_info_string("gene_syn".as_bytes(), &[dbrow.gene_syn.as_bytes()]).expect("failed to set gene_syn field");
    }
    if dbrow.gene_full != "." {
        record.push_info_string("gene_ful".as_bytes(), &[dbrow.gene_full.as_bytes()]).expect("failed to set gene_full field");
    }
    if dbrow.pathway_uniprot != "." {
        record.push_info_string("pathway_uniprot".as_bytes(), &[dbrow.pathway_uniprot.as_bytes()]).expect("failed to set pathway_uniprot field");
    }
    if dbrow.pathway_biocarta != "." {
        record.push_info_string("pathway_biocarta".as_bytes(), &[dbrow.pathway_biocarta.as_bytes()]).expect("failed to set pathway_biocarta field");
    }
    if dbrow.pathway_consensusPathDB != "." {
        record.push_info_string("pathway_consensusPathDB".as_bytes(), &[dbrow.pathway_consensusPathDB.as_bytes()]).expect("failed to set pathway_consensusPathDB field");
    }
    if dbrow.pathway_kegg != "." {
        record.push_info_string("pathway_kegg".as_bytes(), &[dbrow.pathway_kegg.as_bytes()]).expect("failed to set pathway_kegg field");
    }
    if dbrow.gene_function != "." {
        record.push_info_string("gene_function".as_bytes(), &[dbrow.gene_function.as_bytes()]).expect("failed to set gene_function field");
    }
    if dbrow.gene_disease != "." {
        record.push_info_string("gene_disease".as_bytes(), &[dbrow.gene_disease.as_bytes()]).expect("failed to set gene_disease field");
    }
    if dbrow.MIM_phenotype_id != "." {
        record.push_info_string("MIM_phenotype_id".as_bytes(), &[dbrow.MIM_phenotype_id.as_bytes()]).expect("failed to set MIM_phenotype_id field");
    }
    if dbrow.MIM_disease != "." {
        record.push_info_string("MIM_disease".as_bytes(), &[dbrow.MIM_disease.as_bytes()]).expect("failed to set MIM_disease field");
    }
    if dbrow.orphanet_id != "." {
        record.push_info_string("orphanet_id".as_bytes(), &[dbrow.orphanet_id.as_bytes()]).expect("failed to set orphanet_id field");
    }
    if dbrow.orphanet_disorder != "." {
        record.push_info_string("orphanet_disorder".as_bytes(), &[dbrow.orphanet_disorder.as_bytes()]).expect("failed to set orphanet_disorder field");
    }
    if dbrow.orphanet_assoc_type != "." {
        record.push_info_string("orphanet_assoc_type".as_bytes(), &[dbrow.orphanet_assoc_type.as_bytes()]).expect("failed to set orphanet_assoc_type field");
    }
    if dbrow.GWAS_trait != "." {
        record.push_info_string("GWAS_trait".as_bytes(), &[dbrow.GWAS_trait.as_bytes()]).expect("failed to set GWAS_trait field");
    }
    if dbrow.HPO_id != "." {
        record.push_info_string("HPO_id".as_bytes(), &[dbrow.HPO_id.as_bytes()]).expect("failed to set HPO_id field");
    }
    if dbrow.HPO_name != "." {
        record.push_info_string("HPO_name".as_bytes(), &[dbrow.HPO_name.as_bytes()]).expect("failed to set HPO_name field");
    }
    if dbrow.GO_bio_process != "." {
        record.push_info_string("GO_bio_process".as_bytes(), &[dbrow.GO_bio_process.as_bytes()]).expect("failed to set GO_bio_process field");
    }
    if dbrow.GO_cellular_comp != "." {
        record.push_info_string("GO_cellular_comp".as_bytes(), &[dbrow.GO_cellular_comp.as_bytes()]).expect("failed to set GO_cellular_comp field");
    }
    if dbrow.GO_molecular_func != "." {
        record.push_info_string("GO_molecular_func".as_bytes(), &[dbrow.GO_molecular_func.as_bytes()]).expect("failed to set GO_molecular_func field");
    }
    if dbrow.UNIPROT_tissue_specificity != "." {
        record.push_info_string("UNIPROT_tissue_specificity".as_bytes(), &[dbrow.UNIPROT_tissue_specificity.as_bytes()]).expect("failed to set UNIPROT_tissue_specificity field");
    }
    if dbrow.egenetics_expression != "." {
        record.push_info_string("egenetics_expression".as_bytes(), &[dbrow.egenetics_expression.as_bytes()]).expect("failed to set egenetics_expression field");
    }
    if dbrow.GNF_atlas_expression != "." {
        record.push_info_string("GNF_atlas_expression".as_bytes(), &[dbrow.GNF_atlas_expression.as_bytes()]).expect("failed to set GNF_atlas_expression field");
    }
    if dbrow.MGI_mouse_gene != "." {
        record.push_info_string("MGI_mouse_gene".as_bytes(), &[dbrow.MGI_mouse_gene.as_bytes()]).expect("failed to set MGI_mouse_gene field");
    }
    if dbrow.MGI_mouse_phenotype != "." {
        record.push_info_string("MGI_mouse_phenotype".as_bytes(), &[dbrow.MGI_mouse_phenotype.as_bytes()]).expect("failed to set MGI_mouse_phenotype field");
    }
    return
}

pub fn glkup(input: Option<&str>, output: Option<&str>, fields: Option<&str>, dbnsfp: &bool, threads: &usize) {
    let mut bcf = bcfutils::get_rdr(input);
    bcf.set_threads(threads.clone()).expect("unable to set reader threads");

    let hdrv = bcf.header();
    let mut hdr = Header::from_template(&hdrv);

    hdr.push_record(format!("##INFO=<ID=gnomad_gene,Number=1,Type=String,Description=\"gene used to get gnomad gene level fields, pulled from {} INFO fields in descending priority\">", fields.expect("must supply fields to pull gene from")).as_bytes());
    hdr.push_record(r#"##INFO=<ID=gnomAD_pLI,Number=1,Type=Float,Description="gnomAD pLI for gene in gnomad_gene INFO field">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=oe_lof,Number=1,Type=Float,Description="gnomAD oe_lof for gene in gnomad_gene INFO field">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=oe_lof_upper,Number=1,Type=Float,Description="gnomAD oe_lof_upper for gene in gnomad_gene INFO field">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=syn_z,Number=1,Type=Float,Description="gnomAD syn_z for gene in gnomad_gene INFO field">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=mis_z,Number=1,Type=Float,Description="gnomAD mis_z for gene in gnomad_gene INFO field">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=lof_z,Number=1,Type=Float,Description="gnomAD lof_z for gene in gnomad_gene INFO field">"#.as_bytes());
    hdr.push_record(r#"##INFO=<ID=exac_pLI,Number=1,Type=Float,Description="gnomAD exac_pLI for gene in gnomad_gene INFO field">"#.as_bytes());

    let mut obcf = bcfutils::get_wrtr(output, &hdr);
    obcf.set_threads(threads.clone()).expect("unable to set writer threads");

    let fs: Vec<String> = match fields {
        Some(f) => f.split(",").map(|s| s.to_string()).collect(),
        None => {
            eprintln!("Error: need to specify INFO fields to grab gene names from");
            process::exit(1);
        }
    }; 

    let gmap = build_gnomad_map();
    let dmap = build_dbnsfp_map();

    let mut b = Buffer::new();
    for record_result in bcf.records() {
        let mut record = record_result.expect("fail to read record");
        obcf.translate(&mut record);

        for f in &fs {
            let gene = match record.info_shared_buffer(f.as_bytes(), &mut b).string().unwrap() {
                Some(g) => str::from_utf8(g[0]).unwrap(),
                None => continue,
            };
            let grow = match gmap.get(gene) {
                Some(g) => g,
                None => continue,
            };
            add_gnomad_fields(grow, &mut record);
        };

        for f in &fs {
            let gene = match record.info_shared_buffer(f.as_bytes(), &mut b).string().unwrap() {
                Some(g) => str::from_utf8(g[0]).unwrap(),
                None => continue,
            };
            let dbrow = match dmap.get(gene) {
                Some(d) => d,
                None => continue,
            };
            add_dbnsfp_fields(dbrow, &mut record);
        };



        obcf.write(&record).expect("failed to write record");
    }
}
