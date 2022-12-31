#[cfg(target_os="linux")]
pub mod cspice {
    // use std::ffi::CString;
    use std::os::raw::c_char;
    use std::os::raw::c_double;

    #[link(name = "cspice")]
    extern  {
        #[allow(dead_code)]
        // #[allow(improper_ctypes)]
        fn spkpos_c(targ: *const c_char, et: c_double, _ref: *const c_char, abcorr: *const c_char, obs: *const c_char, ptarg: *mut [c_double; 3], it: *mut c_double);

        fn furnsh_c(file: *const c_char);
    }

    #[allow(dead_code)]
    pub fn call_splpos_c_moon_earth_j2000(_t:f64, pos: &mut [f64;3]) {
        // let mut pos:[f64; 3] = [0.0, 0.0, 0.0];
        let mut it:f64 = 0.0;
        unsafe {
            spkpos_c(
                String::from("MOON\0").as_ptr() as *const i8,
                _t,
                String::from("J2000\0").as_ptr() as *const i8,
                String::from("NONE\0").as_ptr() as *const i8,
                String::from("EARTH\0").as_ptr() as *const i8,
                pos,
                &mut it
            );
        }
    }

    #[allow(dead_code)]
    pub fn call_splpos_c_sun_earth_j2000(_t:f64, pos: &mut [f64;3]) {
        // let mut pos:[f64; 3] = [0.0, 0.0, 0.0];
        let mut it:f64 = 0.0;
        unsafe {
            spkpos_c(
                String::from("SUN\0").as_ptr() as *const i8,
                _t,
                String::from("J2000\0").as_ptr() as *const i8,
                String::from("NONE\0").as_ptr() as *const i8,
                String::from("EARTH\0").as_ptr() as *const i8,
                pos,
                &mut it
            );
        }
    }

    #[allow(dead_code)]
    pub fn call_furnsh_c () {
        unsafe {
            furnsh_c(String::from("spk_data/de421.bsp").as_ptr() as *const i8);
        }
    }
}