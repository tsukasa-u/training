#![allow(dead_code, non_camel_case_types, unused_imports, unused_assignments)]

use super::gmp::mp_limb_t;
use super::gmp::mpq_srcptr;
use super::gmp::mpz_srcptr;
use super::mpfr::mpfr_prec_t;
use super::mpfr::mpfr_srcptr;
use super::mpfi::mpfi_srcptr;
use super::mpfr::mpfr_t;
use super::mpfi::__mpfi_struct;
use super::*;

use std::ffi::c_ulong;
use std::ffi::c_char;

pub type mpfi_t = __mpfi_struct;

impl mpfi_t {
    pub fn new() -> Self {
        let mut ret: Self = Self {
            left: mpfr_t::new(),
            right: mpfr_t::new()
        };
        ret.init();
        return ret;
    }

    pub fn new2(prec: u64) -> Self {
        let mut ret: Self = Self {
            left: mpfr_t::new(),
            right: mpfr_t::new()
        };
        ret.init2(prec);
        return ret;
    }

    // Dont't use it !! Use init2
    pub fn init(&mut self) {
        unsafe {
            mpfi::mpfi_init(self);
        }
    }

    pub fn init2(&mut self, prec: u64) {
        unsafe {
            mpfi::mpfi_init2(self, prec as mpfr_prec_t);
        }
    }

    pub fn clear(&mut self) {
        unsafe {
            mpfi::mpfi_clear(self);
        }
    }

    pub fn get_prec(&self) -> u64 {
        let mut ret: mpfr_prec_t = 0;
        unsafe {
            ret = mpfi::mpfi_get_prec(self);
        }
        return ret as u64;
    }

    pub fn set_prec(&mut self, prec: u64) {
        unsafe {
            mpfi::mpfi_set_prec(self, prec as mpfr_prec_t);
        }
    }

    pub fn set_str(&mut self, op: &str) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfi::mpfi_set_str(self, (*op).as_ptr() as *const c_char);
        }
        return ret;
    }

    //  get double value which is returned from function mpfi_mid
    pub fn get_d(&self) -> f64 {
        let mut ret: f64 = 0.0;
        unsafe {
            ret = mpfi::mpfi_get_d(self);
        }
        return ret;
    }

    pub fn get_fr(&self) -> mpfr_t {
        let mut ret: mpfr_t = mpfr_t::new() ;
        unsafe {
            mpfi::mpfi_get_fr(&mut ret, self);
        }
        return ret;
    }

    pub fn revert_if_needed(&mut self) -> i32 {
        let mut ret: i32 = 0;
        unsafe {
            ret = mpfi::mpfi_revert_if_needed(self);
        }
        return ret;
    }

}

macro_rules! fn0_const {
    ($s1: ident $s2: ident) => (
        pub fn $s1() -> mpfi_t {
            let mut tmp: mpfi_t = mpfi_t::new();
            unsafe {
                mpfi::$s2(&mut tmp);
            }
            return tmp;
        }
    )
}

macro_rules! fn1_ex {
    ($s1: ident $s2: ident | $t: ty | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(op1: $t, op2: mpfi_srcptr) -> $u {
                let mut ret: $u = 0 as $u;
                unsafe {
                    ret = mpfi::$s2(op1, op2);
                }
                return ret;
            }
        }
    )
}

macro_rules! fn0_impl {
    ($s1: ident $s2: ident | | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&self) -> $u {
                let mut tmp: $u = <$u>::new();
                unsafe {
                    mpfi::$s2(&mut tmp, self);
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
                let mut tmp: $u = <$u>::new();
                unsafe {
                    mpfi::$s2(&mut tmp, self, x);
                }
                return tmp;
            }
        }
    )
}

macro_rules! fn2_impl {
    ($s1: ident $s2: ident | $t: ty | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&mut self, op: $t) -> $u {
                let mut ret: $u = 0 as $u;
                unsafe {
                    ret = mpfi::$s2(self, op);
                }
                return ret;
            }
        }
    );
    ($s1: ident $s2: ident | $t1: ty, $t2: ty | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&mut self, op1: $t1, op2: $t2) -> $u {
                let mut ret: $u = 0 as $u;
                unsafe {
                    ret = mpfi::$s2(self, op1, op2);
                }
                return ret;
            }
        }
    );
}

macro_rules! fn4_impl {
    ($s1: ident $s2: ident | $t: ty | $u: ty) => (
        impl mpfi_t {
            pub fn $s1(&self, x: $t) -> $u {
                let mut tmp: $u = 0 as $u;
                unsafe {
                    tmp = mpfi::$s2(self, x);
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
                let mut tmp: $u = 0 as $u;
                unsafe {
                    tmp = mpfi::$s2(self);
                }
                return tmp;
            }
        }
    )
}

fn2_impl! {set      mpfi_set    | mpfi_srcptr   | i32   }
fn2_impl! {set_si   mpfi_set_si     | i64   | i32   }
fn2_impl! {set_ui   mpfi_set_ui     | u64   | i32   }
fn2_impl! {set_d    mpfi_set_d      | f64   | i32   }
fn2_impl! {set_flt  mpfi_set_flt    | f32   | i32   }
fn2_impl! {set_fr   mpfi_set_fr     | mpfr_srcptr   | i32   }

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
fn1_impl! { atan2   mpfi_atan2  | mpfi_srcptr   | mpfi_t    }

fn0_impl! { sech    mpfi_sech   |   | mpfi_t    }
fn0_impl! { csch    mpfi_csch   |   | mpfi_t    }
fn0_impl! { coth    mpfi_coth   |   | mpfi_t    }

fn0_impl! { log1p   mpfi_log1p      |   | mpfi_t    }
fn0_impl! { log10p1 mpfi_log10p1    |   | mpfi_t    }
fn0_impl! { log2p1  mpfi_log2p1     |   | mpfi_t    }
fn0_impl! { expm1   mpfi_expm1      |   | mpfi_t    }
fn0_impl! { exp2m1  mpfi_exp2m1     |   | mpfi_t    }
fn0_impl! { exp10m1 mpfi_exp10m1    |   | mpfi_t    }

fn0_impl! { log2    mpfi_log2   |   | mpfi_t    }
fn0_impl! { log10   mpfi_log10  |   | mpfi_t    }

fn1_impl! { hypot   mpfi_hypot  | mpfi_srcptr   | mpfi_t    }

fn4_impl! { cmp     mpfi_cmp        | mpfi_srcptr   | i32   }

fn4_impl! { cmp_d   mpfi_cmp_d      | f64   | i32   }
fn4_impl! { cmp_ui  mpfi_cmp_ui     | u64   | i32   }
fn4_impl! { cmp_si  mpfi_cmp_si     | i64   | i32   }

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

fn0_impl! { get_left    mpfi_get_left   |   | mpfr_t    }
fn0_impl! { get_right   mpfi_get_right  |   | mpfr_t    }

fn2_impl! { put     mpfi_put    | mpfi_srcptr   | i32   }
fn2_impl! { put_si  mpfi_put_si     | i64   | i32   }
fn2_impl! { put_ui  mpfi_put_ui     | u64   | i32   }
fn2_impl! { put_d   mpfi_put_d      | f64   | i32   }
fn2_impl! { put_fr  mpfi_put_fr     | mpfr_srcptr   | i32   }

fn2_impl! { interv_si   mpfi_interv_si  | i64, i64    | i32   }
fn2_impl! { interv_ui   mpfi_interv_ui  | u64, u64    | i32   }
fn2_impl! { interv_d    mpfi_interv_d   | f64, f64    | i32   }
fn2_impl! { interv_fr   mpfi_interv_fr  | mpfr_srcptr, mpfr_srcptr    | i32   }

fn5_impl! { is_empty    mpfi_is_empty   |   | i32   }
fn2_impl! { intersect   mpfi_intersect  | mpfi_srcptr, mpfi_srcptr    | i32   }
fn2_impl! { union       mpfi_union      | mpfi_srcptr, mpfi_srcptr    | i32   }

macro_rules! op2_impl {
    ($s1: ident $s2: ident $s3: ident | $($t: ty)*) => ($(
        impl std::ops::$s1<$t> for mpfi_t {
            type Output =  mpfi_t;

            fn $s2(self, rhs: $t) -> Self::Output {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s3(&mut tmp, &self, rhs);
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

            fn $s2(self, rhs: mpfi_t) -> Self::Output {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s3(&mut tmp, self, &rhs);
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

            fn $s2(self, rhs: mpfi_t) -> Self::Output {
                let mut tmp: mpfi_t = mpfi_t::new();
                unsafe {
                    mpfi::$s3(&mut tmp, &rhs, self);
                }
                return tmp;
            }
        }
    )*)
}

op2_impl! { Add add mpfi_add | mpfi_srcptr }
op2_impl! { Sub sub mpfi_sub | mpfi_srcptr }
op2_impl! { Mul mul mpfi_mul | mpfi_srcptr }
op2_impl! { Div div mpfi_div | mpfi_srcptr }

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
            mpfi::mpfi_neg(&mut tmp, &self);
        }
        return tmp;
    }
}


pub fn swap(x: &mut mpfi_t, y: &mut mpfi_t) {
    unsafe {
        mpfi::mpfi_swap(x, y);
    }
}

fn0_const! { const_log2      mpfi_const_log2     }
fn0_const! { const_pi        mpfi_const_pi       }
fn0_const! { const_euler     mpfi_const_euler    }
fn0_const! { const_catalan   mpfi_const_catalan  }

fn1_ex! { is_strictly_inside mpfi_is_strictly_inside | mpfi_srcptr   | i32   }
fn1_ex! { is_inside          mpfi_is_inside          | mpfi_srcptr   | i32   }
fn1_ex! { is_inside_d        mpfi_is_inside_d        | f64           | i32   }
fn1_ex! { is_inside_ui       mpfi_is_inside_ui       | u64           | i32   }
fn1_ex! { is_inside_si       mpfi_is_inside_si       | i64           | i32   }
fn1_ex! { is_inside_fr       mpfi_is_inside_fr       | mpfr_srcptr   | i32   }

pub fn error() -> i32 {
    let mut tmp: i32 = 0;
    unsafe {
        tmp = mpfi::mpfi_error();
    }
    return tmp;
}

pub fn reset_error() {
    unsafe {
        mpfi::mpfi_reset_error();
    }
}

pub fn set_error(x: i32) {
    unsafe {
        mpfi::mpfi_set_error(x);
    }
}

pub fn is_error() -> i32 {
    let mut tmp: i32 = 0;
    unsafe {
        tmp = mpfi::mpfi_is_error();
    }
    return tmp;
}