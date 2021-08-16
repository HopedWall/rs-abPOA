use ab_poa::abpoa::*;
use std::ptr;
use std::os::raw::c_int;

fn main() {
    unsafe {
        let ab = abpoa_init();
        let mut abpt = abpoa_init_para();

        (*abpt).set_out_msa(1);
        (*abpt).set_out_cons(1);
        (*abpt).w = 6;
        (*abpt).k = 9;
        (*abpt).min_w = 10;
        (*abpt).set_progressive_poa(1);

        abpoa_post_set_para(abpt);

        let seqs: Vec<&str> = [
            "CGTCAATCTATCGAAGCATACGCGGGCAGAGCCGAAGACCTCGGCAATCCA",
            "CCACGTCAATCTATCGAAGCATACGCGGCAGCCGAACTCGACCTCGGCAATCAC",
            "CGTCAATCTATCGAAGCATACGCGGCAGAGCCCGGAAGACCTCGGCAATCAC",
            "CGTCAATGCTAGTCGAAGCAGCTGCGGCAGAGCCGAAGACCTCGGCAATCAC",
            "CGTCAATCTATCGAAGCATTCTACGCGGCAGAGCCGACCTCGGCAATCAC",
            "CGTCAATCTAGAAGCATACGCGGCAAGAGCCGAAGACCTCGGCCAATCAC",
            "CGTCAATCTATCGGTAAAGCATACGCTCTGTAGCCGAAGACCTCGGCAATCAC",
            "CGTCAATCTATCTTCAAGCATACGCGGCAGAGCCGAAGACCTCGGCAATC",
            "CGTCAATGGATCGAGTACGCGGCAGAGCCGAAGACCTCGGCAATCAC",
            "CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT"
        ].to_vec();

        let _nt4_table : [u8; 256] = [
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
            4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
        ];

        // Get the number of input sequences
        let n_seqs: c_int = seqs.len() as c_int;

        // Create a Vec with the sequences' length
        let mut seq_lens: Vec<c_int> = seqs.iter().map(|s| s.len() as c_int).collect();
        //let seq_lens: *mut c_int  = seq_lens_val.as_mut_ptr();
        //println!("seq lens: {:#?}", seq_lens_val);

        // Create a matrix (bseqs_val) where:
        // - each row represents an input sequence (rows can have different lengths)
        // - each column represents a base as per the _nth4_table
        // Bseqs is a Vec of pointers each pointing to one row

        let mut bseqs_val_val : Vec<Vec<u8>> = seqs.into_iter().map(|s| {
            s.chars()
                .map(|c| *(_nt4_table).get(c as usize).unwrap())
                .collect()
        }).collect();
        //println!("bseqs val val: {:#?}",bseqs_val_val);

        let mut bseqs_val: Vec<*mut u8> = bseqs_val_val.iter_mut().map(|s| s.as_mut_ptr()).collect();
        //let bseqs : *mut *mut u8 = bseqs_val.as_mut_ptr();

        // Now perform the alignment
        let mut cons_seq: *mut *mut u8 = ptr::null_mut();
        let mut cons_c: *mut *mut c_int = ptr::null_mut();
        let mut cons_l : *mut c_int = ptr::null_mut();
        let mut cons_n: c_int = 0;
        let mut msa_seq: *mut *mut u8 = ptr::null_mut();
        let mut msa_l : c_int = 0;
        let out : *mut FILE = ptr::null_mut(); //stdout;

        abpoa_msa(ab, abpt, n_seqs, ptr::null_mut(), seq_lens.as_mut_ptr(), bseqs_val.as_mut_ptr(), out,
                  &mut cons_seq, &mut cons_c, &mut cons_l,
                  &mut cons_n, &mut msa_seq, &mut msa_l);

        let alphabet = String::from("ACGTN-");
        let alphabet_vec : Vec<char> = alphabet.chars().collect();
        println!("msa_l: {}", msa_l);
        for i in 0..n_seqs {
            let outer_pointer = *msa_seq.add((i) as usize);
            for j in 0..msa_l {
                let inner_pointer = *(outer_pointer.add(j as usize));
                print!("{}", alphabet_vec.get(inner_pointer as usize).unwrap());
            }
            println!("");
        }

        //println!("out: {}", (*out));
        //println!("Msa_seq: {:#?}", msa_seq);
    }
}
