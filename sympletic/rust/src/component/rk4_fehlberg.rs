
use crate::component::{Vec3::Vec3, VecL::VecL};

const N: usize =  1 + 1 + 2 ;   // satelite (x + v) + moon x + sun x

#[allow(dead_code)]
fn single_step(
    _t: f64,
    _h: f64,
    func: fn(f64, &VecL<Vec3<f64>, N>) -> VecL<Vec3<f64>, N>,
    _x: &VecL<Vec3<f64>, N>,
    _x_pp: &mut VecL<Vec3<f64>, N>,
    #[allow(non_snake_case)]
    Delta_x: &mut VecL<Vec3<f64>, N>
) {
    let k1:VecL<Vec3<f64>, N> = func(_t, _x);
    let k2:VecL<Vec3<f64>, N> = func(_t + _h/4.0,          &(*_x + k1*_h/4.0));
    let k3:VecL<Vec3<f64>, N> = func(_t + _h*3.0/8.0,      &(*_x + (k1*3.0 + k2*9.0)*_h/32.0));
    let k4:VecL<Vec3<f64>, N> = func(_t + _h*12.0/13.0,    &(*_x + (k1*1932.0 - k2*7200.0 + k3*7296.0)*_h/2197.0));
    let k5:VecL<Vec3<f64>, N> = func(_t + _h,              &(*_x + (k1*439.0/216.0 - k2*8.0 + k3*3680.0/513.0 - k4*845.0/4104.0)*_h));
    let k6:VecL<Vec3<f64>, N> = func(_t + _h/2.0,          &(*_x + (-k1*8.0/27.0 + k2*2.0 - k3*3544.0/2565.0 + k4*1859.0/4104.0 - k5*11.0/40.0)*_h));

    *_x_pp = *_x + (k1*16.0/135.0 + k3*6656.0/12825.0 + k4*28561.0/56430.0 - k5*9.0/50.0 + k6*2.0/55.0)*_h;
    // let _x_pp_:Vec3 = _x + _h*(25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0);
    *Delta_x = *_x + (k1/360.0 - k3*128.0/4275.0 - k4*2197.0/75240.0 + k5/50.0 + k6*2.0/55.0)*_h;
}

#[allow(dead_code)]
pub fn rk4_fehlberg(
    _t: &mut f64,
    _h: f64,
    func: fn(f64, &VecL<Vec3<f64>, N>) -> VecL<Vec3<f64>, N>,
    _x: &mut VecL<Vec3<f64>, N>,
    eps:f64
) {

    let mut _x_pp: VecL<Vec3<f64>, N> = VecL::new();
    #[allow(non_snake_case)]
    let mut Delta_x: VecL<Vec3<f64>, N> = VecL::new();
    single_step(*_t, _h, func, _x, &mut _x_pp, &mut Delta_x);

    let absmax = Delta_x.iter().reduce(|accum, item| {
        if accum.distance2() >= item.distance2() { accum } else { item }
    })
        .unwrap()
        .distance2()
        .powf(0.5);

    let beta = 0.9*(eps/absmax).powf(0.2);

    if beta > 1.0 {
        *_x = _x_pp;
    } else {
        let ndiv: f64 = if 10.0 - 1.0/beta <= 0.0 { 10.0 } else { (1.0/beta).floor() };
        
        for i in 0..ndiv as usize {
            single_step(*_t + (i as f64)/ndiv, _h/ndiv, func, _x, &mut _x_pp, &mut Delta_x);
            *_x = _x_pp;
        }
    }
    *_t += _h
}

// [satelite_x, satelite_v, moon_a, sun_x]
fn funcVecL(_t: f64, _x:&VecL<Vec3<f64>, N>) -> VecL<Vec3<f64>, N> {
    return VecL::new();
}