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

            // c_int is an alias for i32
            let n_seqs: c_int = 10;
            //let seq_lens: *mut c_int  = malloc((mem::size_of::<i32>() * n_seqs as usize) as u64) as *mut c_int;
            let mut seq_lens_val : [c_int; 10] = [0; 10];

            let seq_lens: *mut c_int  = seq_lens_val.as_mut_ptr();

            //let bseqs: *mut *mut u8 = malloc((mem::size_of::<u8>() * n_seqs as usize) as u64) as *mut *mut u8;
            let mut bseqs_val: [* mut u8; 10] = [ptr::null_mut();10];
            let bseqs = bseqs_val.as_mut_ptr();
            for i in 0..n_seqs {
                let curr_seq : &str = seqs.get(i as usize).unwrap();
                //seq_lens[i] = &strlen(seqs[i]);
                let mut curr_seq_len : &mut c_int = seq_lens_val.get_mut(i as usize).unwrap();
                *curr_seq_len = curr_seq.len() as c_int;
                //bseqs[i] = malloc(u8::size_of() * seq_lens[i]);
                //println!("Curr seq: {:#?}", curr_seq);
                //let curr_seq_as_chars = curr_seq.chars();
                //println!("Chars: {:#?}\n", curr_seq_as_chars);
                let mut curr_numbers : Vec<u8> = curr_seq.chars()
                                                .map(|c| *(_nt4_table).get(c as usize).unwrap())
                                                .collect();
                //println!("Curr numbers: {:#?}", curr_numbers);
                /*
                for char in curr_seq_as_chars {
                    println!("Char {} is as usize: {} value matrix: {} (as u8: {})", char, char as usize, *(_nt4_table).get(char as usize).unwrap(), *(_nt4_table).get(char as usize).unwrap() as u8);
                }
                 */

                let mut curr_bseq_val = bseqs_val.get_mut(i as usize).unwrap();
                *curr_bseq_val = curr_numbers.as_mut_ptr();
            }

            /*
            println!("ab: {:#?}", ab);
            println!("abpt: {:#?}", abpt);
            println!("Seq lens: {:#?}", seq_lens_val);
            for i in 0..n_seqs {
                let mut curr_bseq_val = bseqs_val.get_mut(i as usize).unwrap();

                let curr_seq_lens = *seq_lens_val.get(i as usize).unwrap();
                println!("Curr seq lens: {}", curr_seq_lens);
                for j in 0..curr_seq_lens {
                    println!("Bseq [{}][{}]: {:#?}", i, j, *(curr_bseq_val.add(j as usize)));
                }

                println!("");
            }


            abpoa_msa(ab, abpt, n_seqs, ptr::null_mut(), seq_lens, bseqs, stdout,
                      ptr::null_mut(), ptr::null_mut(), ptr::null_mut(),
                      ptr::null_mut(), ptr::null_mut(), ptr::null_mut());
            */

            /*
            let cons_seq: *mut *mut *mut u8 = ptr::null_mut();
            let cons_cov: *mut *mut *mut c_int = ptr::null_mut();
            let cons_l : *mut *mut c_int = ptr::null_mut();
            let cons_n: *mut c_int = &mut 0;
            let msa_seq: *mut *mut *mut u8 = ptr::null_mut();
            let msa_l : *mut c_int = &mut 0;

            abpoa_msa(ab_poa, abpt, 10, ptr::null_mut(), seq_lens, bseqs, ptr::null_mut(), cons_seq, cons_cov, cons_l, cons_n, msa_seq, msa_l);

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
