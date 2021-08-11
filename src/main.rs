mod abpoa;

use abpoa::*;

fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::mem;

    #[test]
    fn test() {
        unsafe {
            let ab_poa = abpoa_init();
            let mut abpt = abpoa_init_para();

            /*
            //(*abpt).out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
            abpt.as_ref.unwrap().set_out_msa(1);
            //(*abpt).out_cons = 1; // generate consensus sequence, set 0 to disable
            abpt.as_ref_mut().unwrap().set_out_cons(1);
            (*abpt).w = 6;
            (*abpt).k = 9;
            (*abpt).min_w = 10; // minimizer-based seeding and partition
            //(*abpt).progressive_poa = 1;
            abpt.as_ref_mut().unwrap().set_progressive_poa(1);
             */

            abpoa_post_set_para(abpt);
        }
    }
}
