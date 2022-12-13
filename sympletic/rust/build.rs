use std::env;
// use std::path::PathBuf;
extern crate bindgen;

fn main(){
    let project_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    println!("cargo:rustc-link-search={}/cspice/lib/", project_dir);
    println!("cargo:rustc-link-lib=static=cspice");

}