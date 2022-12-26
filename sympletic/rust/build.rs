use std::env;
// use std::path::PathBuf;
extern crate bindgen;

fn main(){
    let project_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    println!("cargo:rustc-link-search={}/cspice/lib/", project_dir);
    println!("cargo:rustc-link-lib=static=cspice");

    // println!("cargo:rustc-link-search={}/mpfi-1.5.4/src/.libs", project_dir);
    // println!("cargo:rustc-link-lib=init");
    
    // println!("cargo:rustc-link-search={}/mpfi-1.5.4/src", project_dir);
    // println!("cargo:rustc-link-lib=init");
    
    println!("cargo:rustc-link-search=/usr/local/lib");
    println!("cargo:rustc-link-lib=mpfi")
    
}