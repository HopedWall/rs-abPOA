# rs-abPOA
Rust bindings for the [abPOA partial order aligner](https://github.com/yangao07/abPOA). 
Only works under *NIX systems.

## Contents
There are two modules in this crate:

- **abpoa.rs** contains the raw bindings, generated with [bindgen](https://github.com/rust-lang/rust-bindgen)
- **abpoa_wrapper.rs** contains a "Rust-friendly" wrapper to interact with abPOA

## Usage
### Alignment
```
let mut aligner = AbpoaAligner::new_with_example_params();

let seqs: Vec<&str> = [
    "CGTCAAT",
    "CCACGTCAAT",
    "CGTCAAT",
    "CGTCAATGCTA",
].to_vec();

let aln = aligner.align_seqs(&seqs);
```

### Consensus
```
let mut aligner = AbpoaAligner::new_with_example_params();

let seqs: Vec<&str> = [
    "CGTCAAT",
    "CCACGTCAAT",
    "CGTCAAT",
    "CGTCAATGCTA",
].to_vec();

let cons = aligner.consensus_from_seqs(&seqs);
```

## Interacting with abPOA
This crate also provides an abstraction layer over abPOA, 
that allows for the definition of nodes labelled with seqs of any size. 
```
let mut aligner = AbpoaAligner::new_with_example_params();

// Create this graph: "ACG" -> "AAA"
aligner.add_nodes_edges(&vec!["ACG", "AAA"], &vec![(0, 1)]);

// Align "ACG" to the graph defined above
let res = aligner.align_sequence("ACG");
```

## TODOs
- [X] Raw bindings and wrapper
- [X] Find a more efficient way to abpoa_id -> graph_id
- [ ] Build abPOA from source
- [ ] Limit unsafe usage
