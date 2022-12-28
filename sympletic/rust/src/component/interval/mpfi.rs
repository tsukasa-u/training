#![allow(dead_code, non_camel_case_types, unused_imports)]

use super::gmp::gmp_randstate_t;
use super::gmp::mpq_srcptr;
use super::gmp::mpz_srcptr;
use super::mpfr::mpfr_t;
use super::mpfr::mpfr_ptr;
use super::mpfr::mpfr_srcptr;

use std::ffi::c_void;
use std::ffi::c_long;
use std::ffi::c_int;
use std::ffi::c_ulong;
use std::ffi::c_double;
use std::ffi::c_float;
use std::ffi::c_char;
use std::ffi::c_longlong;
use std::ffi::c_ulonglong;

#[repr(C)]
pub struct __mpfi_struct {
    pub left: mpfr_t,
    pub right: mpfr_t
}

pub type mpfi_ptr = *mut __mpfi_struct;
pub type mpfi_srcptr = *const __mpfi_struct;
type mpfr_prec_t = c_long;

#[link(name="mpfi")]
extern {
    pub fn mpfi_round_prec  (x: mpfi_ptr,  prec: mpfr_prec_t    )   -> c_void;

    pub fn mpfi_init        (x: mpfi_ptr                        )   -> c_void;
    pub fn mpfi_init2       (x: mpfi_ptr,   prec: mpfr_prec_t   )   -> c_void;
    pub fn mpfi_inits       (x: mpfi_ptr,   ...                 )   -> c_void;
    pub fn mpfi_inits2      (p: mpfr_prec_t, x: mpfi_ptr, ...   )   -> c_void;

    pub fn mpfi_clear       (x: mpfi_ptr        )   -> c_void;
    pub fn mpfi_clears      (x: mpfi_ptr, ...   )   -> c_void;
    
    pub fn mpfi_get_prec    (x: mpfi_srcptr                 )   -> mpfr_prec_t;
    pub fn mpfi_set_prec    (x: mpfi_ptr, prec: mpfr_prec_t )   -> c_void;

    /* assignment functions                         */
    pub fn mpfi_set         (a: mpfi_ptr, b: mpfi_srcptr    )   -> c_int;
    pub fn mpfi_set_si      (a: mpfi_ptr, b: c_long         )   -> c_int;
    pub fn mpfi_set_ui      (a: mpfi_ptr, b: c_ulong        )   -> c_int;
    pub fn mpfi_set_d       (a: mpfi_ptr, b: c_double       )   -> c_int;
    pub fn mpfi_set_flt     (a: mpfi_ptr, b: c_float        )   -> c_int;
    // int     mpfi_set_ld     (mpfi_ptr, const long double);
    pub fn mpfi_set_z       (a: mpfi_ptr, b: mpz_srcptr     )   -> c_int;
    pub fn mpfi_set_q       (a: mpfi_ptr, b: mpq_srcptr     )   -> c_int;
    pub fn mpfi_set_fr      (a: mpfi_ptr, b: mpfr_srcptr    )   -> c_int;
    pub fn mpfi_set_str     (a: mpfi_ptr, b: *const c_char  )   -> c_int;

    /* int     mpfi_set_sj     (mpfi_ptr, const intmax_t); */
    /* int     mpfi_set_uj     (mpfi_ptr, const uintmax_t); */

    /* combined initialization and assignment functions */
    pub fn mpfi_init_set    (a: mpfi_ptr, b: mpfi_srcptr    )   -> c_int;
    pub fn mpfi_init_set_si (a: mpfi_ptr, b: c_long         )   -> c_int;
    pub fn mpfi_init_set_ui (a: mpfi_ptr, b: c_ulong        )   -> c_int;
    pub fn mpfi_init_set_d  (a: mpfi_ptr, b: c_double       )   -> c_int;
    pub fn mpfi_init_set_flt(a: mpfi_ptr, b: c_float        )   -> c_int;
    // int     mpfi_init_set_ld    (mpfi_ptr, const long double);
    pub fn mpfi_init_set_z  (a: mpfi_ptr, b: mpz_srcptr     )   -> c_int;
    pub fn mpfi_init_set_q  (a: mpfi_ptr, b: mpq_srcptr     )   -> c_int;
    pub fn mpfi_init_set_fr (a: mpfi_ptr, b: mpfr_srcptr    )   -> c_int;
    pub fn mpfi_init_set_str(a: mpfi_ptr, b: *const c_char  )   -> c_int;

    /* swapping two intervals */
    pub fn mpfi_swap        (a: mpfi_ptr, b: mpfi_ptr   )   -> c_void;


    /* Various useful interval functions            */
    /* with scalar or interval results              */

    /* absolute diameter                            */
    pub fn mpfi_diam_abs    (diam: mpfr_ptr, interv: mpfi_srcptr)   -> c_int;
    /* relative diameter                            */
    pub fn mpfi_diam_rel    (diam: mpfr_ptr, interv: mpfi_srcptr)   -> c_int;
    /* diameter: relative if the interval does not contain 0 */
    /* absolute otherwise                                    */
    pub fn mpfi_diam        (diam: mpfr_ptr, interv: mpfi_srcptr)   -> c_int;
    /* magnitude: the largest absolute value of any element */
    pub fn mpfi_mag         (m: mpfr_ptr, y: mpfi_srcptr)   -> c_int;
    /* mignitude: the smallest absolute value of any element */
    pub fn mpfi_mig         (m: mpfr_ptr, y: mpfi_srcptr)   -> c_int;
    /* middle of y                                           */
    pub fn mpfi_mid         (m: mpfr_ptr, y: mpfi_srcptr)   -> c_int;
    /* picks randomly a point m in y */
    pub fn mpfi_alea        (m: mpfr_ptr, y: mpfi_srcptr)   -> c_void;

    pub fn mpfi_urandom     (m: mpfr_ptr, y: mpfi_srcptr, state: gmp_randstate_t)   -> c_void;
    pub fn mpfi_nrandom     (m: mpfr_ptr, y: mpfi_srcptr, state: gmp_randstate_t)   -> c_void;
    pub fn mpfi_erandom     (m: mpfr_ptr, y: mpfi_srcptr, state: gmp_randstate_t)   -> c_void;


    /* Conversions                                  */

    pub fn mpfi_get_d       (a: mpfi_srcptr             )   -> c_double;
    pub fn mpfi_get_fr      (m: mpfr_ptr, a: mpfi_srcptr)   -> c_void;


    /* Basic arithmetic operations                  */

    /* arithmetic operations between two interval operands */
    pub fn mpfi_add         (a: mpfi_ptr, b: mpfi_srcptr, c: mpfi_srcptr)   -> c_int;
    pub fn mpfi_sub         (a: mpfi_ptr, b: mpfi_srcptr, c: mpfi_srcptr)   -> c_int;
    pub fn mpfi_mul         (a: mpfi_ptr, b: mpfi_srcptr, c: mpfi_srcptr)   -> c_int;
    pub fn mpfi_div         (a: mpfi_ptr, b: mpfi_srcptr, c: mpfi_srcptr)   -> c_int;

    /* arithmetic operations between an interval operand and a double prec. floating-point */
    pub fn mpfi_add_d       (a: mpfi_ptr, b: mpfi_srcptr, c: c_double    )   -> c_int;
    pub fn mpfi_sub_d       (a: mpfi_ptr, b: mpfi_srcptr, c: c_double    )   -> c_int;
    pub fn mpfi_d_sub       (a: mpfi_ptr, b: c_double   , c: mpfi_srcptr )   -> c_int;
    pub fn mpfi_mul_d       (a: mpfi_ptr, b: mpfi_srcptr, c: c_double    )   -> c_int;
    pub fn mpfi_div_d       (a: mpfi_ptr, b: mpfi_srcptr, c: c_double    )   -> c_int;
    pub fn mpfi_d_div       (a: mpfi_ptr, b: c_double   , c: mpfi_srcptr )   -> c_int;

    /* arithmetic operations between an interval operand and an unsigned long integer */
    pub fn mpfi_add_ui      (a: mpfi_ptr, b: mpfi_srcptr, c: c_ulong     )   -> c_int;
    pub fn mpfi_sub_ui      (a: mpfi_ptr, b: mpfi_srcptr, c: c_ulong     )   -> c_int;
    pub fn mpfi_ui_sub      (a: mpfi_ptr, b: c_ulong    , c: mpfi_srcptr )   -> c_int;
    pub fn mpfi_mul_ui      (a: mpfi_ptr, b: mpfi_srcptr, c: c_ulong     )   -> c_int;
    pub fn mpfi_div_ui      (a: mpfi_ptr, b: mpfi_srcptr, c: c_ulong     )   -> c_int;
    pub fn mpfi_ui_div      (a: mpfi_ptr, b: c_ulong    , c: mpfi_srcptr )   -> c_int;

    // /* arithmetic operations between an interval operand and a long integer */
    pub fn mpfi_add_si      (a: mpfi_ptr, b: mpfi_srcptr, c: c_long      )   -> c_int;
    pub fn mpfi_sub_si      (a: mpfi_ptr, b: mpfi_srcptr, c: c_long      )   -> c_int;
    pub fn mpfi_si_sub      (a: mpfi_ptr, b: c_long     , c: mpfi_srcptr )   -> c_int;
    pub fn mpfi_mul_si      (a: mpfi_ptr, b: mpfi_srcptr, c: c_long      )   -> c_int;
    pub fn mpfi_div_si      (a: mpfi_ptr, b: mpfi_srcptr, c: c_long      )   -> c_int;
    pub fn mpfi_si_div      (a: mpfi_ptr, b: c_long     , c: mpfi_srcptr )   -> c_int;

    /* arithmetic operations between an interval operand and a multiple prec. integer */
    pub fn mpfi_add_z       (a: mpfi_ptr, b: mpfi_srcptr, c: mpz_srcptr  )   -> c_int;
    pub fn mpfi_sub_z       (a: mpfi_ptr, b: mpfi_srcptr, c: mpz_srcptr  )   -> c_int;
    pub fn mpfi_z_sub       (a: mpfi_ptr, b: mpz_srcptr , c: mpfi_srcptr )   -> c_int;
    pub fn mpfi_mul_z       (a: mpfi_ptr, b: mpfi_srcptr, c: mpz_srcptr  )   -> c_int;
    pub fn mpfi_div_z       (a: mpfi_ptr, b: mpfi_srcptr, c: mpz_srcptr  )   -> c_int;
    pub fn mpfi_z_div       (a: mpfi_ptr, b: mpz_srcptr , c: mpfi_srcptr )   -> c_int;

    /* arithmetic operations between an interval operand and a multiple prec. rational */
    pub fn mpfi_add_q       (a: mpfi_ptr, b: mpfi_srcptr, c: mpq_srcptr  )   -> c_int;
    pub fn mpfi_sub_q       (a: mpfi_ptr, b: mpfi_srcptr, c: mpq_srcptr  )   -> c_int;
    pub fn mpfi_q_sub       (a: mpfi_ptr, b: mpq_srcptr , c: mpfi_srcptr )   -> c_int;
    pub fn mpfi_mul_q       (a: mpfi_ptr, b: mpfi_srcptr, c: mpq_srcptr  )   -> c_int;
    pub fn mpfi_div_q       (a: mpfi_ptr, b: mpfi_srcptr, c: mpq_srcptr  )   -> c_int;
    pub fn mpfi_q_div       (a: mpfi_ptr, b: mpq_srcptr , c: mpfi_srcptr )   -> c_int;

    /* arithmetic operations between an interval operand and a mult. prec. floating-pt nb */
    pub fn mpfi_add_fr      (a: mpfi_ptr, b: mpfi_srcptr, c: mpfr_srcptr )   -> c_int;
    pub fn mpfi_sub_fr      (a: mpfi_ptr, b: mpfi_srcptr, c: mpfr_srcptr )   -> c_int;
    pub fn mpfi_fr_sub      (a: mpfi_ptr, b: mpfr_srcptr, c: mpfi_srcptr )   -> c_int;
    pub fn mpfi_mul_fr      (a: mpfi_ptr, b: mpfi_srcptr, c: mpfr_srcptr )   -> c_int;
    pub fn mpfi_div_fr      (a: mpfi_ptr, b: mpfi_srcptr, c: mpfr_srcptr )   -> c_int;
    pub fn mpfi_fr_div      (a: mpfi_ptr, b: mpfr_srcptr, c: mpfi_srcptr )   -> c_int;

    /* arithmetic operations taking a single interval operand */
    pub fn mpfi_neg         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_sqr         (a: mpfi_ptr, u: mpfi_srcptr )   -> c_int;   //  square ?
    /* the inv function generates the whole real interval
    if 0 is in the interval defining the divisor */
    pub fn mpfi_inv         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;   // inverse ?
    /* the sqrt of a (partially) negative interval is a NaN */
    pub fn mpfi_sqrt        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_cbrt        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    /* the first interval contains the absolute values of */
    /* every element of the second interval */
    pub fn mpfi_abs         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    // int     mpfi_inv        (mpfi_ptr, mpfi_srcptr);
    pub fn mpfi_rec_sqrt    (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;

    /* extended division: returns 2 intervals if the denominator contains 0 */
    pub fn mpfi_div_ext	    (res1: mpfi_ptr, res2: mpfi_ptr, op1: mpfi_srcptr, op2: mpfi_srcptr)   -> c_int;

    /* various operations */
    pub fn mpfi_mul_2exp    (a: mpfi_ptr, b: mpfi_srcptr, c: c_ulong )   -> c_int;
    pub fn mpfi_mul_2ui     (a: mpfi_ptr, b: mpfi_srcptr, c: c_ulong )   -> c_int;
    pub fn mpfi_mul_2si     (a: mpfi_ptr, b: mpfi_srcptr, c: c_long  )   -> c_int;
    pub fn mpfi_div_2exp    (a: mpfi_ptr, b: mpfi_srcptr, c: c_ulong )   -> c_int;
    pub fn mpfi_div_2ui     (a: mpfi_ptr, b: mpfi_srcptr, c: c_ulong )   -> c_int;
    pub fn mpfi_div_2si     (a: mpfi_ptr, b: mpfi_srcptr, c: c_long  )   -> c_int;

    // /* Special functions                                        */
    pub fn mpfi_log         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_exp         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_exp2        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_exp10       (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;

    pub fn mpfi_cos         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_sin         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_tan         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_acos        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_asin        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_atan        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_atan2       (a: mpfi_ptr, b: mpfi_srcptr, c: mpfi_srcptr )   -> c_int;

    pub fn mpfi_sec         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_csc         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_cot         (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;

    pub fn mpfi_cosh        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_sinh        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_tanh        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_acosh       (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_asinh       (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_atanh       (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;

    pub fn mpfi_sech        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_csch        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_coth        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;

    pub fn mpfi_log1p       (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_log10p1     (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_log2p1      (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_expm1       (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_exp2m1      (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_exp10m1     (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;

    pub fn mpfi_log2        (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;
    pub fn mpfi_log10       (a: mpfi_ptr, b: mpfi_srcptr )   -> c_int;

    pub fn mpfi_hypot       (a: mpfi_ptr, b: mpfi_srcptr, c: mpfi_srcptr )   -> c_int;

    pub fn mpfi_const_log2      (a: mpfi_ptr)  -> c_int;
    pub fn mpfi_const_pi        (a: mpfi_ptr)  -> c_int;
    pub fn mpfi_const_euler     (a: mpfi_ptr)  -> c_int;
    pub fn mpfi_const_catalan   (a: mpfi_ptr)  -> c_int;

    // /* Comparison functions                                     */
    // /* Warning: the meaning of interval comparison is not clearly defined */
    // /* customizable comparison functions */

    // extern int    (*mpfi_cmp)       (mpfi_srcptr, mpfi_srcptr);

    // extern int    (*mpfi_cmp_d)     (mpfi_srcptr, const double);
    // extern int    (*mpfi_cmp_ui)    (mpfi_srcptr, const unsigned long);
    // extern int    (*mpfi_cmp_si)    (mpfi_srcptr, const long);
    // extern int    (*mpfi_cmp_z)     (mpfi_srcptr, mpz_srcptr);
    // extern int    (*mpfi_cmp_q)     (mpfi_srcptr, mpq_srcptr);
    // extern int    (*mpfi_cmp_fr)    (mpfi_srcptr, mpfr_srcptr);

    // extern int    (*mpfi_is_pos)    (mpfi_srcptr);
    // extern int    (*mpfi_is_nonneg) (mpfi_srcptr);
    // extern int    (*mpfi_is_neg)    (mpfi_srcptr);
    // extern int    (*mpfi_is_nonpos) (mpfi_srcptr);
    // extern int    (*mpfi_is_zero)   (mpfi_srcptr);
    // extern int    (*mpfi_is_strictly_pos) (mpfi_srcptr);
    // extern int    (*mpfi_is_strictly_neg) (mpfi_srcptr);


    pub fn mpfi_cmp         (a: mpfi_srcptr, b: mpfi_srcptr )   -> c_int;

    pub fn mpfi_cmp_d       (a: mpfi_srcptr, b: c_double    )   -> c_int;
    pub fn mpfi_cmp_ui      (a: mpfi_srcptr, b: c_ulong     )   -> c_int;
    pub fn mpfi_cmp_si      (a: mpfi_srcptr, b: c_long      )   -> c_int;
    pub fn mpfi_cmp_z       (a: mpfi_srcptr, b: mpz_srcptr  )   -> c_int;
    pub fn mpfi_cmp_q       (a: mpfi_srcptr, b: mpq_srcptr  )   -> c_int;
    pub fn mpfi_cmp_fr      (a: mpfi_srcptr, b: mpfr_srcptr )   -> c_int;

    pub fn mpfi_is_pos      (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_nonneg   (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_neg      (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_nonpos   (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_zero     (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_strictly_pos     (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_strictly_neg     (a: mpfi_srcptr)   -> c_int;
    

    pub fn mpfi_has_zero    (a: mpfi_srcptr)   -> c_int;

    pub fn mpfi_nan_p       (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_inf_p       (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_bounded_p   (a: mpfi_srcptr)   -> c_int;

    /* Interval manipulation */

    /* operations related to the internal representation by endpoints */

    /* get left or right bound of the interval defined by the
    second argument and put the result in the first one */
    pub fn mpfi_get_left    (b: mpfr_ptr, a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_get_right   (b: mpfr_ptr, a: mpfi_srcptr)   -> c_int;

    pub fn mpfi_revert_if_needed    (a: mpfi_ptr)   -> c_int;

    /* Set operations on intervals */
    /* "Convex hulls" */
    /* extends the interval defined by the first argument
    so that it contains the second one */

    pub fn mpfi_put         (a: mpfi_ptr, b: mpfi_srcptr)   -> c_int;
    pub fn mpfi_put_d       (a: mpfi_ptr, b: c_double   )   -> c_int;
    pub fn mpfi_put_si      (a: mpfi_ptr, b: c_long     )   -> c_int;
    pub fn mpfi_put_ui      (a: mpfi_ptr, b: c_ulong    )   -> c_int;
    pub fn mpfi_put_z       (a: mpfi_ptr, b: mpz_srcptr )   -> c_int;
    pub fn mpfi_put_q       (a: mpfi_ptr, b: mpq_srcptr )   -> c_int;
    pub fn mpfi_put_fr      (a: mpfi_ptr, b: mpfr_srcptr)   -> c_int;

    /* builds an interval whose left bound is the lower (round -infty)
    than the second argument and the right bound is greater
    (round +infty) than the third one */

    pub fn mpfi_interv_d    (a: mpfi_ptr, b: c_double   , c: c_double    )   -> c_int;
    pub fn mpfi_interv_si   (a: mpfi_ptr, b: c_long     , c: c_long      )   -> c_int;
    pub fn mpfi_interv_ui   (a: mpfi_ptr, b: c_ulong    , c: c_ulong     )   -> c_int;
    pub fn mpfi_interv_z    (a: mpfi_ptr, b: mpz_srcptr , c: mpz_srcptr  )   -> c_int;
    pub fn mpfi_interv_q    (a: mpfi_ptr, b: mpq_srcptr , c: mpq_srcptr  )   -> c_int;
    pub fn mpfi_interv_fr   (a: mpfi_ptr, b: mpfr_srcptr, c: mpfr_srcptr )   -> c_int;

    /* Inclusion tests */
    /* tests if the first argument is inside the interval
    defined by the second one */
    pub fn mpfi_is_strictly_inside  (a: mpfi_srcptr, b: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_inside        	(a: mpfi_srcptr, b: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_inside_d      	(a: c_double   , b: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_inside_ui     	(a: c_ulong    , b: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_inside_si     	(a: c_long     , b: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_inside_z      	(a: mpz_srcptr , b: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_inside_q      	(a: mpq_srcptr , b: mpfi_srcptr)   -> c_int;
    pub fn mpfi_is_inside_fr   	    (a: mpfr_srcptr, b: mpfi_srcptr)   -> c_int;

    /* set operations */
    pub fn mpfi_is_empty            (a: mpfi_srcptr)   -> c_int;
    pub fn mpfi_intersect           (a: mpfi_ptr, b: mpfi_srcptr, c: mpfi_srcptr)   -> c_int;
    pub fn mpfi_union               (a: mpfi_ptr, b: mpfi_srcptr, c: mpfi_srcptr)   -> c_int;

    /* complement... : to do later */


    /* Miscellaneous */

    /* adds the second argument to the right bound of the first one
    and subtracts the second argument to the left bound of
    the first one */
    pub fn mpfi_increase    (a: mpfi_ptr, e: mpfr_srcptr)   -> c_int;
    /* keeps the same center and multiply the radius by 2*(1+fact) */
    pub fn mpfi_blow        (y: mpfi_ptr, x: mpfi_srcptr, fact: c_double)   -> c_int;
    /* splits the interval into 2 halves */
    pub fn mpfi_bisect      (y1: mpfi_ptr, y2: mpfi_ptr, y: mpfi_srcptr)   -> c_int;

    pub fn mpfi_get_version (           )   -> *const c_char;

    /* Error handling */

    pub fn mpfi_error       (           )   -> c_int;
    pub fn mpfi_reset_error (           )   -> c_void;
    pub fn mpfi_set_error   (i: c_int   )   -> c_void;
    pub fn mpfi_is_error    (           )   -> c_int;
}