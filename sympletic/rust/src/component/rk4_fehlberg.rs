
use crate::component::{Vec3::Vec3, VecL::VecL};

const N: usize = 2;

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
    let k1:VecL<Vec3<f64>, N> = func(_t, &(*_x).clone());
    let k2:VecL<Vec3<f64>, N> = func(_t + _h/4.0,          &((*_x).clone() + k1.clone()*_h/4.0));
    let k3:VecL<Vec3<f64>, N> = func(_t + _h*3.0/8.0,      &((*_x).clone() + (k1.clone()*3.0 + k2.clone()*9.0)*_h/32.0));
    let k4:VecL<Vec3<f64>, N> = func(_t + _h*12.0/13.0,    &((*_x).clone() + (k1.clone()*1932.0 - k2.clone()*7200.0 + k3.clone()*7296.0)*_h/2197.0));
    let k5:VecL<Vec3<f64>, N> = func(_t + _h,              &((*_x).clone() + (k1.clone()*439.0/216.0 - k2.clone()*8.0 + k3.clone()*3680.0/513.0 - k4.clone()*845.0/4104.0)*_h));
    let k6:VecL<Vec3<f64>, N> = func(_t + _h/2.0,          &((*_x).clone() + (-k1.clone()*8.0/27.0 + k2.clone()*2.0 - k3.clone()*3544.0/2565.0 + k4.clone()*1859.0/4104.0 - k5.clone()*11.0/40.0)*_h));

    *_x_pp = (*_x).clone() + (k1.clone()*16.0/135.0 + k3.clone()*6656.0/12825.0 + k4.clone()*28561.0/56430.0 - k5.clone()*9.0/50.0 + k6.clone()*2.0/55.0)*_h;
    // let _x_pp_:Vec3 = _x + _h*(25.0*k1/216.0 + 1408.0*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0);
    *Delta_x = (*_x).clone() + (k1/360.0 - k3*128.0/4275.0 - k4*2197.0/75240.0 + k5/50.0 + k6*2.0/55.0)*_h;
}

#[allow(dead_code)]
pub fn rk4_fehlberg() {

}