
#![allow(dead_code, non_camel_case_types, unused_imports, unused_assignments)]

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
use std::ffi::CStr;


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

pub type mpfr_t = __mpfr_struct;

impl mpfr_t {
    pub fn new() -> Self {
        return Self {
            _mpfr_prec : 0,
            _mpfr_sign : 0,
            _mpfr_exp :  0,
            _mpfr_d :  0 as *mut mp_limb_t
        };
    }
}

#[repr(C)]
pub enum mpfr_rnd_t {
    MPFR_RNDN=0,  /* round to nearest, with ties to even */
    MPFR_RNDZ,    /* round toward zero */
    MPFR_RNDU,    /* round toward +Inf */
    MPFR_RNDD,    /* round toward -Inf */
    MPFR_RNDA,    /* round away from zero */
    MPFR_RNDF,    /* faithful rounding */
    MPFR_RNDNA=-1 /* round to nearest, with ties away from zero (mpfr_round) */
}

/* Compatibility with previous types of MPFR */


// type mp_rnd_t = mpfr_rnd_t;
// type mp_prec_t = mpfr_prec_t;

pub type mpfr_ptr = *mut __mpfr_struct;
pub type mpfr_srcptr = *const __mpfr_struct;

#[link(name="mpfr")]
extern {

    pub fn mpfr_get_exp(x: mpfr_srcptr) -> mpfr_exp_t;
    pub fn mpfr_set_exp(x: mpfr_ptr, y: mpfr_exp_t) -> c_int;
    pub fn mpfr_get_prec(x: mpfr_srcptr) -> mpfr_prec_t;
    pub fn mpfr_set_prec(x: mpfr_ptr, y: mpfr_prec_t) -> c_void;
    pub fn mpfr_set_default_prec(x: mpfr_prec_t) -> c_void;
    pub fn mpfr_get_default_prec() -> mpfr_prec_t;

    pub fn mpfr_get_flt(op: mpfr_srcptr, rnd: mpfr_rnd_t) -> c_float;
    pub fn mpfr_get_d(op: mpfr_srcptr, rnd: mpfr_rnd_t) -> c_double;
    pub fn mpfr_get_d1(op: mpfr_srcptr) -> c_double;
    pub fn mpfr_get_d_2exp(y: *mut c_long, op: mpfr_srcptr, rnd: mpfr_rnd_t) -> c_double;
    pub fn mpfr_frexp(y: *mut mpfr_exp_t, op: mpfr_srcptr, rad: mpfr_rnd_t) -> c_long;
    pub fn mpfr_get_si(op: mpfr_srcptr, rad: mpfr_rnd_t) -> c_long;
    pub fn mpfr_get_ui(op: mpfr_srcptr, rad: mpfr_rnd_t) -> c_ulong;
    pub fn mpfr_get_str(
        str: *mut c_char,           // The result string.
        expptr: *mut mpfr_exp_t,    // The returned exponent.
        b: c_int,                   // The base which vary from 2 to 62. 
        n: usize,                   // The number of digits in the result string.
        op: mpfr_srcptr,
        rnd: mpfr_rnd_t
    ) -> *const c_char;             // Return the converted string of digits.
}

impl mpfr_t {
    pub fn mpfr_get_exp(&self) -> mpfr_exp_t {
        let mut ret: mpfr_exp_t = 0;
        unsafe {
            ret = mpfr_get_exp(self);
        }
        return ret;
    }

    pub fn mpfr_set_exp(&mut self, a: mpfr_exp_t) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfr_set_exp(self, a);
        }
        return ret;
    }
    
    pub fn mpfr_get_prec(&self) -> mpfr_prec_t {
        let mut ret: mpfr_prec_t = 0;
        unsafe {
            ret = mpfr_get_prec(self);
        }
        return ret;
    }
    
    pub fn mpfr_set_prec(&mut self, a: mpfr_prec_t) {
        unsafe {
            mpfr_set_prec(self, a);
        }
    }
    
    pub fn mpfr_get_flt(&self, rnd: mpfr_rnd_t) -> c_float {
        let mut ret: c_float = 0.0;
        unsafe {
            ret = mpfr_get_flt(self, rnd);
        }
        return ret;
    }
    
    pub fn mpfr_get_d(&self, rnd: mpfr_rnd_t) -> c_double {
        let mut ret: c_double = 0.0;
        unsafe {
            ret = mpfr_get_d(self, rnd);
        }
        return ret;
    }

    pub fn mpfr_get_d1(&self) -> c_double {
        let mut ret: c_double = 0.0;
        unsafe {
            ret = mpfr_get_d1(self);
        }
        return ret;
    }
    
    pub fn mpfr_get_d_2exp(&self, y: *mut c_long, rnd: mpfr_rnd_t) -> c_double {
        let mut ret: c_double = 0.0;
        unsafe {
            ret = mpfr_get_d_2exp(y,self,  rnd);
        }
        return ret;
    }
    
    pub fn mpfr_frexp(&self, y: *mut mpfr_exp_t, rnd: mpfr_rnd_t) -> c_long {
        let mut ret: c_long = 0;
        unsafe {
            ret = mpfr_frexp(y,self,  rnd);
        }
        return ret;
    }
    
    pub fn mpfr_get_si(&self, rnd: mpfr_rnd_t) -> c_long {
        let mut ret: c_long = 0;
        unsafe {
            ret = mpfr_get_si(self, rnd);
        }
        return ret;
    }
    
    pub fn mpfr_get_ui(&self, rnd: mpfr_rnd_t) -> c_ulong {
        let mut ret: c_ulong = 0;
        unsafe {
            ret = mpfr_get_ui(self, rnd);
        }
        return ret;
    }

    pub fn get_str(&self, expptr: &mut mpfr_exp_t, b: c_int, n: usize, rnd: mpfr_rnd_t) -> &str {
        let mut ret: &str = "";
        let str_ptr: *mut c_char = 0 as *mut i8;
        unsafe {
            let tmp: &CStr = CStr::from_ptr(
                mpfr_get_str(str_ptr, expptr, b, n, self, rnd)
            );
            ret = tmp.to_str().unwrap();
        }
        return ret;
    }

    pub fn get_str2(&self, expptr: &mut mpfr_exp_t, b: c_int, n: usize, rnd: mpfr_rnd_t) -> &str {
        let mut ret: &str = "";
        let str_ptr: *mut c_char = 0 as *mut i8;
        unsafe {
            mpfr_get_str(str_ptr, expptr, b, n, self, rnd);
            ret = CStr::from_ptr(str_ptr).to_str().unwrap();
        }
        return ret;
    }

}

pub fn set_default_prec(x: mpfr_prec_t) {
    unsafe {
        mpfr_set_default_prec(x);
    }
}

pub fn get_default_prec() -> mpfr_prec_t {
    let mut ret: mpfr_prec_t = 0;
    unsafe {
        ret = mpfr_get_default_prec();
    }
    return ret;
}