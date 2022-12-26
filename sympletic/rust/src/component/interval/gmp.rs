
#![allow(dead_code, non_camel_case_types)]


#[allow(unused_imports)]
use std::ffi::c_longlong;
#[allow(unused_imports)]
use std::ffi::c_ulonglong;
use std::ffi::c_long;
use std::ffi::c_ulong;
use std::ffi::c_int;
use std::ffi::c_void;

#[cfg(target_pointer_width="32")]
type mp_limb_t = c_uint;
#[cfg(target_pointer_width="32")]
type mp_limb_signed_t = c_int;

// #[cfg(target_pointer_width="64")]
//  type mp_limb_t = c_ulonglong;
// #[cfg(target_pointer_width="64")]
//  type mp_limb_signed_t = c_longlong;
#[cfg(target_pointer_width="64")]
pub type mp_limb_t = c_ulong;
#[cfg(target_pointer_width="64")]
type mp_limb_signed_t = c_long;
    
#[repr(C)]
pub struct __mpz_struct {
    _mp_alloc: c_int,
    _mp_size: c_int ,
    _mp_d: *mut mp_limb_t
    // _mp_d: *mut c_void
}

    
type MP_INT = __mpz_struct ;
type mpz_t = *mut __mpz_struct ;
    
type mp_ptr = *mut mp_limb_t;
type mp_srcptr = *const mp_limb_t;

// TODO: What is _CRAY and _CRAYMPP
//  #if defined (_CRAY) && ! defined (_CRAYMPP)
//  #define __GMP_MP_SIZE_T_INT     1
#[cfg(target_pointer_width="32")]
type mp_size_t = c_int;
#[cfg(target_pointer_width="32")]
type mp_exp_t = c_int;
//  #else
//  #define __GMP_MP_SIZE_T_INT     0
#[cfg(target_pointer_width="64")]
type mp_size_t = c_long;
#[cfg(target_pointer_width="64")]
type mp_exp_t = c_long;
//  #endif
    
#[repr(C)]
pub struct __mpq_struct {
    _mp_num: __mpz_struct,
    _mp_den: __mpz_struct
}
    
type MP_RAT = __mpq_struct;
type mpq_t = *mut __mpq_struct;
    
#[repr(C)]
struct __mpf_struct {
    _mp_prec : c_int,
    _mp_size: c_int,
    _mp_exp: mp_exp_t,
    _mp_d: *mut mp_limb_t
}

type mpf_t = *mut __mpf_struct;
/* Available random number generation algorithms.  */
#[repr(C)]
enum gmp_randalg_t {
    GMP_RAND_ALG_DEFAULT = 0,
    // GMP_RAND_ALG_LC = GMP_RAND_ALG_DEFAULT /* Linear congruential.  */
    GMP_RAND_ALG_LC
} 

#[repr(C)]
union __func_ptr {
    _mp_lc : *mut c_void         /* Pointer to function pointers structure.  */
}

/* Random state struct.  */
#[repr(C)]
pub struct __gmp_randstate_struct {
    _mp_seed: mpz_t,         /* _mp_d member points to state of the generator. */
    _mp_alg: gmp_randalg_t,  /* Currently unused. */
    _mp_algdata: __func_ptr
}

pub type gmp_randstate_t = *mut __gmp_randstate_struct;

/* Types for function declarations in gmp files.  */
pub type mpz_srcptr = *const __mpz_struct;
type mpz_ptr = *mut __mpz_struct;
type mpf_srcptr = *const __mpf_struct;
type mpf_ptr = *mut __mpf_struct;
pub type mpq_srcptr = *const __mpq_struct;
type mpq_ptr = *mut __mpq_struct;