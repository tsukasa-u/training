#![allow(dead_code, non_camel_case_types, unused_imports)]

use super::gmp::mp_limb_t;
use super::mpfr::__mpfr_struct;
use super::mpfi::__mpfi_struct;
use super::*;

use std::ffi::c_ulong;

pub struct Interval(__mpfi_struct);

impl Interval {
    pub fn new() -> Self {
        return Self(
            __mpfi_struct {
                left: __mpfr_struct {
                    _mpfr_prec : 0,
                    _mpfr_sign : 0,
                    _mpfr_exp :  0,
                    _mpfr_d :  0 as *mut mp_limb_t
                },
                right: __mpfr_struct {
                    _mpfr_prec : 0,
                    _mpfr_sign : 0,
                    _mpfr_exp :  0,
                    _mpfr_d :  0 as *mut mp_limb_t
                }
            }
        );
    }

    pub fn init(&mut self) {
        unsafe {
            mpfi::mpfi_init(&mut self.0);
        }
    }
}