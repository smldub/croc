use std::env;
use std::path::PathBuf;
use cmake::Config;
use std::io::Write;

fn main() -> std::io::Result<()>{
    // Build libcint with cmake
    let dst = Config::new("libcint").build();

    // Tell cargo to look for shared libraries in the specified directory
    // You might need to change it from lib64 to something else
    println!("cargo:rustc-link-search={}/lib64", dst.display());

    // Tell cargo to tell rustc to link the system 
    // shared library.
    println!("cargo:rustc-link-lib=cint");

    // Need to overwrite a header file because libcint assumes
    // that we installed the library (when we didn't hehe)
    let contents = std::fs::read_to_string(format!("{}/include/cint_funcs.h",dst.display()))?;
    let new = contents.replace(r##"#include <cint.h>"##,r##"#include "cint.h""##);
    let mut file = std::fs::OpenOptions::new().write(true).truncate(true).open(format!("{}/include/cint_funcs.h",dst.display()))?;
    file.write(new.as_bytes())?;

    // The bindgen::Builder is the main entry point
    // to bindgen, and lets you build up options for
    // the resulting bindings.
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header(format!("{}/include/cint_funcs.h", dst.display()))
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    // Write the bindings to the $OUT_DIR/bindings.rs file.
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    dbg!(bindings)
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
    Ok(())
}
