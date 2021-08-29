extern crate bindgen;

use std::path::PathBuf;

fn main() {
    //Example command to compile .c
    //gcc -g example.c -I ./include -L ./lib -labpoa -lz -o example

    // include abpoa (-L ./lib -labpoa)
    println!("cargo:rustc-link-search=./abpoa-lib");
    println!("cargo:rustc-link-lib=abpoa");

    // include z library (-lz)
    println!("cargo:rustc-link-lib=z");

    // TODO (write in readme): it may be necessary to install clang if not present
    //sudo apt-get install -y clang

    // Tell cargo to invalidate the built crate whenever the wrapper changes
    //println!("cargo:rerun-if-changed=wrapper.h");

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header("wrapper.h")
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from("./src/");
    bindings
        .write_to_file(out_path.join("abpoa.rs"))
        .expect("Couldn't write bindings!");
}
