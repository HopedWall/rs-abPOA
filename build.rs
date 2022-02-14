extern crate bindgen;

use std::path::PathBuf;
use std::process::Command;

fn main() {

    // use make to build abPOA
    Command::new("make")
        .arg("-C")
        .arg("abPOA")
        .output()
        .expect("failed to invoke make");

    // include abpoa (-L ./lib -labpoa)
    println!("cargo:rustc-link-search=abPOA/lib");
    println!("cargo:rustc-link-lib=abpoa");

    // include z library (-lz)
    println!("cargo:rustc-link-lib=z");

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
    let out_path = PathBuf::from("src/");
    bindings
        .write_to_file(out_path.join("abpoa.rs"))
        .expect("Couldn't write bindings!");
}
