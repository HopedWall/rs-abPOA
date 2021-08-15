use ab_poa::abpoa::*;
use std::ptr::null;

fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{mem, io, ptr};
    use std::os::raw::{c_int, c_char};
    use std::alloc::alloc;
    use std::ops::Deref;
    use core::slice;
    use std::mem::MaybeUninit;

    #[test]
    fn test() {
        unsafe {
            let ab = abpoa_init();
            let mut abpt = abpoa_init_para();

            (*abpt).set_out_msa(1);
            (*abpt).set_out_cons(1);
            (*abpt).w = 6;
            (*abpt).k = 9;
            (*abpt).min_w = 10;
            (*abpt).set_progressive_poa(1);

            //println!("abpt: {:#?}", *abpt);

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
            let mut seq_lens_val : Vec<c_int> = seqs.iter().map(|s| s.len() as c_int).collect();
            let seq_lens: *mut c_int  = seq_lens_val.as_mut_ptr();
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

            let mut bseqs_val: Vec<*mut u8> = bseqs_val_val.iter_mut().map(|mut s| s.as_mut_ptr()).collect();
            let bseqs : *mut *mut u8 = bseqs_val.as_mut_ptr();

            // Now perform the alignment
            let cons_seq: *mut *mut *mut u8 = ptr::null_mut();
            let cons_c: *mut *mut *mut c_int = ptr::null_mut();
            let cons_l : *mut *mut c_int = ptr::null_mut();
            let cons_n: *mut c_int = ptr::null_mut();
            let msa_seq: *mut *mut *mut u8 = ptr::null_mut();
            let msa_l : *mut c_int = ptr::null_mut();
            let out : *mut FILE = ptr::null_mut();

            abpoa_msa(ab, abpt, n_seqs, ptr::null_mut(), seq_lens, bseqs, ptr::null_mut(),
                      cons_seq, cons_c, cons_l,
                      cons_n, msa_seq, msa_l);

            println!(">Multiple_sequence_alignment\n");
            println!("MSA len is: {:#?}", msa_l);
            println!("Cons n is: {:#?}", cons_n);

            /*
            for i in 0..n_seqs {
                for j in 0..*msa_l {
                    print!("[i:{}, j:{}]: ", i, j);
                    println!("{:#?}", **msa_seq.add((i*n_seqs+j) as usize));
                }
                println!();
            }
             */

            //let v = slice::from_raw_parts(msa_seq, *msa_l as usize).to_vec();
            //println!("{:#?}",msa_seq);

            //let result : *mut *mut u8 = *msa_seq;
            //let result2 = &result;
            /*
            println!(">Multiple_sequence_alignment");
            for i in 0..n_seqs {
                print!("ACGTN-");
                for j in 0..*msa_l {
                    print!("{}", msa_seq[i][j])
                }
                println!();
            }
             */
        }
    }
}
