
// use gmp_mpfr_sys::mpfr::mpfr_t;
use crate::component::interval::interval;

pub fn tester() {
    let mut x: interval::mpfi_t = interval::mpfi_t::new();
    x.init();
    x.set_si(3);
    let y: interval::mpfi_t = 1 as i64 + x;
    println!("{}", y.get_d())
}