#![allow(dead_code, non_camel_case_types, unused_imports, unused_assignments)]

use super::gmp::mp_limb_t;
use super::gmp::mpq_srcptr;
use super::gmp::mpz_srcptr;
use super::mpfr::mpfr_prec_t;
use super::mpfr::mpfr_srcptr;
use super::mpfi::mpfi_srcptr;
use super::mpfr::__mpfr_struct;
use super::mpfi::__mpfi_struct;
use super::*;

use std::ffi::c_ulong;
use std::ffi::c_char;

pub struct mpfi_t(__mpfi_struct);

impl mpfi_t {
    // Default 53 bits
    pub fn new() -> Self {
        let mut ret: Self = Self(
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
        ret.init();
        return ret;
    }

    pub fn new2(prec: u64) -> Self {
        let mut ret: Self = Self(
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
        ret.init2(prec);
        return ret;
    }

    // Dont't use it !! Use init2
    pub fn init(&mut self) {
        unsafe {
            mpfi::mpfi_init(&mut self.0);
        }
    }

    pub fn init2(&mut self, prec: u64) {
        unsafe {
            mpfi::mpfi_init2(&mut self.0, prec as mpfr_prec_t);
        }
    }

    pub fn clear(&mut self) {
        unsafe {
            mpfi::mpfi_clear(&mut self.0);
        }
    }

    pub fn get_prec(&self) -> u64 {
        let mut ret: mpfr_prec_t = 0;
        unsafe {
            ret = mpfi::mpfi_get_prec(&(*self).0);
        }
        return ret as u64;
    }

    pub fn set_prec(&mut self, prec: u64) {
        unsafe {
            mpfi::mpfi_set_prec(&mut self.0, prec as mpfr_prec_t);
        }
    }

    pub fn set(&mut self, op: &Self) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfi::mpfi_set(&mut self.0, &(*op).0);
        }
        return ret;
    }

    pub fn set_si(&mut self, op: i64) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfi::mpfi_set_si(&mut self.0, op);
        }
        return ret;
    }

    pub fn set_ui(&mut self, op: u64) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfi::mpfi_set_ui(&mut self.0, op);
        }
        return ret;
    }
    
    pub fn set_d(&mut self, op: f64) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfi::mpfi_set_d(&mut self.0, op);
        }
        return ret;
    }
    
    pub fn set_flt(&mut self, op: f32) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfi::mpfi_set_flt(&mut self.0, op);
        }
        return ret;
    }
    
    // pub fn set_fr(&mut self, op: mpfr_srcptr) -> i32 {
    //     let mut ret: i32 = 0;
    //     unsafe {
    //         ret = mpfi::mpfi_set_fr(&mut self.0, op);
    //     }
    //     return ret;
    // }

    pub fn set_str(&mut self, op: &str) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfi::mpfi_set_str(&mut self.0, (*op).as_ptr() as *const c_char);
        }
        return ret;
    }

    //  get double value which is returned from function mpfi_mid
    pub fn get_d(&self) -> f64 {
        let mut ret: f64 = 0.0;
        unsafe {
            ret = mpfi::mpfi_get_d(&self.0);
        }
        return ret;
    }

    // pub fn get_fr(&self) -> Self {
    //     let mut ret: __mpfr_struct ;
    //     unsafe {
    //         mpfi::mpfi_get_fr();
    //     }
    //     return ret;
    // }

}

macro_rules! fn0_const {
    ($s1: ident $s2: ident) => (
        pub fn $s1() -> mpfi_t {
            let mut tmp: mpfi_t = mpfi_t::new();
            unsafe {
                mpfi::$s2(&mut tmp.0);
            }
            return tmp;
        }
    )
}

macro_rules! fn0_impl {
    ($s1: ident $s2: ident | | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&self) -> $u {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s2(&mut tmp.0, &self.0);
                }
                return tmp;
            }
        }
    )
}

macro_rules! fn1_impl {
    ($s1: ident $s2: ident | $t: ty | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&self, x: $t) -> $u {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s2(&mut tmp.0, &self.0, x);
                }
                return tmp;
            }
        }
    )
}


macro_rules! fn3_impl {
    ($s1: ident $s2: ident | $t: ty | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&self, x: & $t) -> $u {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s2(&mut tmp.0, &self.0, &x.0);
                }
                return tmp;
            }
        }
    )
}

macro_rules! fn4_impl {
    ($s1: ident $s2: ident | $t: ty | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&self, x: $t) -> $u {
                let mut tmp: i32 = 0;
                unsafe {
                    tmp = mpfi::$s2(&self.0, x);
                }
                return tmp;
            }
        }
    )
}

macro_rules! fn5_impl {
    ($s1: ident $s2: ident | | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&self) -> $u {
                let mut tmp: i32 = 0;
                unsafe {
                    tmp = mpfi::$s2(&self.0);
                }
                return tmp;
            }
        }
    )
}

fn0_impl! { sqr     mpfi_sqr    |   | mpfi_t    }
fn0_impl! { inv     mpfi_inv    |   | mpfi_t    }
fn0_impl! { sqrt    mpfi_sqrt   |   | mpfi_t    }
fn0_impl! { cbrt    mpfi_cbrt   |   | mpfi_t    }
fn0_impl! { abs     mpfi_abs    |   | mpfi_t    }

fn0_impl! { rec_sqrt    mpfi_rec_sqrt   |   | mpfi_t    }

fn1_impl! { mul_2exp    mpfi_mul_2exp   | u64   | mpfi_t    }
fn1_impl! { mul_2ui     mpfi_mul_2ui    | u64   | mpfi_t    }
fn1_impl! { mul_2si     mpfi_mul_2si    | i64   | mpfi_t    }
fn1_impl! { div_2exp    mpfi_div_2exp   | u64   | mpfi_t    }
fn1_impl! { div_2ui     mpfi_div_2ui    | u64   | mpfi_t    }
fn1_impl! { div_2si    mpfi_div_2si     | i64   | mpfi_t    }

fn0_impl! { log     mpfi_log    |   | mpfi_t    }
fn0_impl! { exp     mpfi_exp    |   | mpfi_t    }
fn0_impl! { exp2    mpfi_exp2   |   | mpfi_t    }
fn0_impl! { exp10   mpfi_exp10  |   | mpfi_t    }
fn0_impl! { cos     mpfi_cos    |   | mpfi_t    }
fn0_impl! { sin     mpfi_sin    |   | mpfi_t    }
fn0_impl! { tan     mpfi_tan    |   | mpfi_t    }
fn0_impl! { acos    mpfi_acos   |   | mpfi_t    }
fn0_impl! { asin    mpfi_asin   |   | mpfi_t    }
fn0_impl! { atan    mpfi_atan   |   | mpfi_t    }
fn3_impl! { atan2   mpfi_atan2  | mpfi_t    | mpfi_t    }

fn0_impl! { sech    mpfi_sech   |   | mpfi_t    }
fn0_impl! { csch    mpfi_csch   |   | mpfi_t    }
fn0_impl! { coth    mpfi_coth   |   | mpfi_t    }

fn0_impl! { log1p    mpfi_log1p     |   | mpfi_t    }
fn0_impl! { log10p1  mpfi_log10p1   |   | mpfi_t    }
fn0_impl! { log2p1   mpfi_log2p1    |   | mpfi_t    }
fn0_impl! { expm1    mpfi_expm1     |   | mpfi_t    }
fn0_impl! { exp2m1   mpfi_exp2m1    |   | mpfi_t    }
fn0_impl! { exp10m1  mpfi_exp10m1   |   | mpfi_t    }

fn0_impl! { log2     mpfi_log2  |   | mpfi_t    }
fn0_impl! { log10    mpfi_log10 |   | mpfi_t    }

fn3_impl! { hypot    mpfi_hypot | mpfi_t    | mpfi_t    }

fn4_impl! { cmp_d    mpfi_cmp_d     | f64   | i32   }
fn4_impl! { cmp_ui   mpfi_cmp_ui    | u64   | i32   }
fn4_impl! { cmp_si   mpfi_cmp_si    | i64   | i32   }

fn5_impl! { pos     mpfi_is_pos     |   | i32   }
fn5_impl! { nonneg  mpfi_is_nonneg  |   | i32   }
fn5_impl! { neg     mpfi_is_neg     |   | i32   }
fn5_impl! { nonpos  mpfi_is_nonpos  |   | i32   }
fn5_impl! { zero    mpfi_is_zero    |   | i32   }
fn5_impl! { strictly_pos    mpfi_is_strictly_pos    |   | i32   }
fn5_impl! { strictly_neg    mpfi_is_strictly_neg    |   | i32   }

fn5_impl! { has_zero    mpfi_has_zero   |   | i32   }


fn5_impl! { nan_p   mpfi_nan_p  |   | i32   }
fn5_impl! { inf_p   mpfi_inf_p  |   | i32   }
fn5_impl! { bounded_p   mpfi_bounded_p  |   | i32   }

macro_rules! op1_impl {
    ($s1: ident $s2: ident $s3: ident | $($t: ty)*) => ($(
        impl std::ops::$s1<$t> for mpfi_t {
            type Output =  mpfi_t;

            // #[inline]
            fn $s2(self, rhs: $t) -> Self::Output {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s3(&mut tmp.0, &self.0, &rhs.0);
                }
                return tmp;
            }
        }
    )*)
}

macro_rules! op2_impl {
    ($s1: ident $s2: ident $s3: ident | $($t: ty)*) => ($(
        impl std::ops::$s1<$t> for mpfi_t {
            type Output =  mpfi_t;

            // #[inline]
            fn $s2(self, rhs: $t) -> Self::Output {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s3(&mut tmp.0, &self.0, rhs);
                }
                return tmp;
            }
        }
    )*)
}

macro_rules! op3_impl {
    ($s1: ident $s2: ident $s3: ident | $($t: ty)*) => ($(
        impl std::ops::$s1<mpfi_t> for $t {
            type Output =  mpfi_t;

            // #[inline]
            fn $s2(self, rhs: mpfi_t) -> Self::Output {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s3(&mut tmp.0, self, &rhs.0);
                }
                return tmp;
            }
        }
    )*)
}

macro_rules! op4_impl {
    ($s1: ident $s2: ident $s3: ident | $($t: ty)*) => ($(
        impl std::ops::$s1<mpfi_t> for $t {
            type Output =  mpfi_t;

            // #[inline]
            fn $s2(self, rhs: mpfi_t) -> Self::Output {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s3(&mut tmp.0, &rhs.0, self);
                }
                return tmp;
            }
        }
    )*)
}

op1_impl! { Add add mpfi_add | mpfi_t }
op1_impl! { Sub sub mpfi_sub | mpfi_t }
op1_impl! { Mul mul mpfi_mul | mpfi_t }
op1_impl! { Div div mpfi_div | mpfi_t }

op2_impl! { Add add mpfi_add_d | f64 }
op4_impl! { Add add mpfi_add_d | f64 }
op2_impl! { Sub sub mpfi_sub_d | f64 }
op3_impl! { Sub sub mpfi_d_sub | f64 }
op2_impl! { Mul mul mpfi_mul_d | f64 }
op4_impl! { Mul mul mpfi_mul_d | f64 }
op2_impl! { Div div mpfi_div_d | f64 }
op3_impl! { Div div mpfi_d_div | f64 }

op2_impl! { Add add mpfi_add_ui | u64 }
op4_impl! { Add add mpfi_add_ui | u64 }
op2_impl! { Sub sub mpfi_sub_ui | u64 }
op3_impl! { Sub sub mpfi_ui_sub | u64 }
op2_impl! { Mul mul mpfi_mul_ui | u64 }
op4_impl! { Mul mul mpfi_mul_ui | u64 }
op2_impl! { Div div mpfi_div_ui | u64 }
op3_impl! { Div div mpfi_ui_div | u64 }

op2_impl! { Add add mpfi_add_si | i64 }
op4_impl! { Add add mpfi_add_si | i64 }
op2_impl! { Sub sub mpfi_sub_si | i64 }
op3_impl! { Sub sub mpfi_si_sub | i64 }
op2_impl! { Mul mul mpfi_mul_si | i64 }
op4_impl! { Mul mul mpfi_mul_si | i64 }
op2_impl! { Div div mpfi_div_si | i64 }
op3_impl! { Div div mpfi_si_div | i64 }

// concat_idents!(mpfi, )

impl std::ops::Neg for mpfi_t {
    type Output = mpfi_t;
    fn neg(self) -> Self::Output {
        let mut tmp: mpfi_t = mpfi_t::new();
        unsafe {
            mpfi::mpfi_neg(&mut tmp.0, &self.0);
        }
        return tmp;
    }
}


pub fn swap(x: &mut mpfi_t, y: &mut mpfi_t) {
    unsafe {
        mpfi::mpfi_swap(&mut (*x).0, &mut (*y).0);
    }
}

fn0_const!{ const_log2      mpfi_const_log2     }
fn0_const!{ const_pi        mpfi_const_pi       }
fn0_const!{ const_euler     mpfi_const_euler    }
fn0_const!{ const_catalan   mpfi_const_catalan  }