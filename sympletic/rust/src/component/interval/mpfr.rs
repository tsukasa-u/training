
#![allow(dead_code, non_camel_case_types, unused_imports)]

use super::gmp::mp_limb_t;

use std::ffi::c_void;
use std::ffi::c_long;
use std::ffi::c_int;
use std::ffi::c_ulong;
use std::ffi::c_double;
use std::ffi::c_float;
use std::ffi::c_char;
use std::ffi::c_longlong;
use std::ffi::c_ulonglong;


#[cfg(target_pointer_width="32")]
type mpfr_prec_t = c_int;
#[cfg(target_pointer_width="32")]
type mpfr_uprec_t = c_uiont;

#[cfg(target_pointer_width="64")]
pub type  mpfr_prec_t = c_long;
#[cfg(target_pointer_width="64")]
pub type mpfr_uprec_t = c_ulong;

/* Definition of sign */
type mpfr_sign_t = c_int;

#[cfg(target_pointer_width="32")]
type mpfr_exp_t = c_int;
#[cfg(target_pointer_width="32")]
type mpfr_uexp_t = c_uint;

#[cfg(target_pointer_width="64")]
type mpfr_exp_t = c_long;
#[cfg(target_pointer_width="64")]
type mpfr_uexp_t = c_ulong;

/* Definition of the main structure */
#[repr(C)]
pub struct __mpfr_struct {
    pub _mpfr_prec: mpfr_prec_t,
    pub _mpfr_sign : mpfr_sign_t,
    pub _mpfr_exp : mpfr_exp_t,
    pub _mpfr_d : *mut mp_limb_t,
}
/* Compatibility with previous types of MPFR */


// type mp_rnd_t = mpfr_rnd_t;
// type mp_prec_t = mpfr_prec_t;

pub type mpfr_ptr = *mut __mpfr_struct;
pub type mpfr_srcptr = *const __mpfr_struct;

// pub struct mpfr_t(__mpfr_struct);