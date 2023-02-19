#![allow(dead_code, unused_imports)]

use rust::test1;
use rust::test2;

// use gmp_mpfr_sys::mpfr::mpfr_t;
use crate::component::interval::interval;
use crate::component::interval::mpfr::mpfr_rnd_t;

macro_rules! stringify_each {
    ($($x:tt)*) => {
        stringify!($($x)/ *)
    }
}

pub fn tester() {
    let x: interval::mpfi_t = interval::mpfi_t::new();
    // // x.init();
    // println!("{:?}", x.get_prec());
    // x.set_prec(120);
    // println!("{:?}", x.get_prec());
    // x.set_si(3);
    // println!("{:?}", &x as *const interval::mpfi_t);
    // println!("{:?}", x);
    // let y: interval::mpfi_t = x;
    // println!("{:?}", &y as *const interval::mpfi_t);
    // println!("{:?}", y);
    // let mut z: interval::mpfi_t = 1 as i64/y;
    // println!("{:?}", z.left.get_d(mpfr_rnd_t::MPFR_RNDN));
    // println!("{:?}", z.right.get_d(mpfr_rnd_t::MPFR_RNDN));
//     // println!("{}", y.get_d());
//     // println!("{}", interval::const_catalan().get_d());
//     // println!("{:?}", interval::const_catalan());
//     // println!("{:?}", interval::const_catalan().get_left());
//     // println!("{:?}", interval::const_catalan().left.mpfr_get_d(mpfr_rnd_t::MPFR_RNDN));
//     // println!("{:?}", interval::const_catalan().right.mpfr_get_d(mpfr_rnd_t::MPFR_RNDN));
    // println!("{:?}", (1 as f64/x.clone()).left.get_d(mpfr_rnd_t::MPFR_RNDN));
    // println!("{:?}", (1 as f64/x.clone()).right.get_d(mpfr_rnd_t::MPFR_RNDN));
    // z.clear();

    // #[test1(test)]
    // fn t() {}

    test2! (
        x in 1..2 {}
    );

    let mut y: interval::mpfi_t = x;

    
    // println!("{}", stringify_each!(abc a b c));
    // println!("{}", stringify_each!(x+x+x+x+x));
    // println!("{}", stringify_each!(&&&&&));
    // println!("{}", stringify_each!(|||||));
    // println!("{}", stringify_each!(<<<<<));
    // println!("{}", stringify_each!(>>>>>));
    // println!("{}", stringify_each!(=====));

    y.clear();
}