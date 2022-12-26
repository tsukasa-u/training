
#[cfg(target_os="linux")]
pub mod interval;

#[cfg(target_os="linux")]
mod gmp;

#[cfg(target_os="linux")]
mod mpfr;

#[cfg(target_os="linux")]
mod mpfi;