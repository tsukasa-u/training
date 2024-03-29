#![allow(dead_code, unused_imports, non_snake_case)]


mod component;

use component::{VecL::Vec3, VecL::VecL, Julian, PhyQty};

mod Moon;

mod c_kernel;
use c_kernel::cspice;

// mod macros;

static INSTANCE_C: once_cell::sync::OnceCell<Vec<f64>> =  once_cell::sync::OnceCell::new();
static INSTANCE_D: once_cell::sync::OnceCell<Vec<f64>> =  once_cell::sync::OnceCell::new();

const PHY_G:f64 = 6.67430E-11;

mod tester;

fn record_orbit(_t: f64, h: f64, interval:f64, satelite: &mut PhyQty::mxvr, _obj: &mut Vec<PhyQty::mxvr>, recoder: &mut Vec<[f32; 10]>) {
    if _t%interval < h {
        let mut ret:[f32;10] = [0.0; 10];
        let mut iter: usize = 4;
        ret[0] = _t as f32;
        ret[1] = (satelite.get_x()[0]/100.0) as f32;
        ret[2] = (satelite.get_x()[1]/100.0) as f32;
        ret[3] = (satelite.get_x()[2]/100.0) as f32;
        for v in _obj {
            if v.get_name().contains("EARTH") == false {
                ret[iter    ] = (v.get_x()[0]/100.0) as f32;
                ret[iter + 1] = (v.get_x()[1]/100.0) as f32;
                ret[iter + 2] = (v.get_x()[2]/100.0) as f32;
                iter += 3;
            }
        }
        (*recoder).push(ret);
    }
}

fn update_orbit(_t: f64, _obj:&mut Vec<PhyQty::mxvr>) {
    for key in _obj {
        match (*key).get_name().as_str() {
            "MOON" => {
                let x = key.get_x_ptr();
                cspice::call_splpos_c_moon_earth_j2000((_t) as f64, (*x).as_mut_ptr());
                *x *= 1000.0;
            },
            "SUN"  => {
                let x = key.get_x_ptr();
                cspice::call_splpos_c_sun_earth_j2000(_t as f64, (*x).as_mut_ptr());
                *x *= 1000.0;
            },
            _      => {}
        }
    }
}

fn Tp(satelite:&PhyQty::mxvr) -> Vec3::Vec3<f64> {
    // println!("{:?}", (*satelite).get_mass());
    return (*satelite).get_v();
}

fn Vq(satelite:&PhyQty::mxvr, obj:&Vec<PhyQty::mxvr>) -> Vec3::Vec3<f64> {
    let mut ret:Vec3::Vec3<f64> = Vec3::Vec3([0.0, 0.0, 0.0]);

    for iter in obj {
        let x_sub: Vec3::Vec3<f64> = (*satelite).get_x() - (*iter).get_x();
        // println!("{} {:?}", iter.get_name(), x_sub/x_sub.distance2());
        ret += x_sub*((*iter).get_mass()/x_sub.distance3());
    }
    // println!("{:?}", ret*PHY_G*satelite.get_mass());
    return ret*PHY_G;
}

fn sympletic4(h:f64, start:f64, end:f64, satelite:&mut PhyQty::mxvr, obj:&mut Vec<PhyQty::mxvr>, recoder: &mut Vec<[f32; 10]>) {

    let c: &Vec<f64> = INSTANCE_C.get().unwrap();
    let d: &Vec<f64> = INSTANCE_C.get().unwrap();

    let mut _t: f64 = start;
    if end - start <= 0.0 {
        return;
    }

    // let x = satelite.get_x_ptr();
    // let v = satelite.get_v_ptr();

    while _t < end {
        update_orbit(_t, obj);
        for iter in 0..4 {
            satelite.x += Tp(&satelite)*(*c)[iter]*h;
            satelite.v -= Vq(&satelite, &obj)*(*d)[iter]*h;
            // println!("{:?}", satelite.get_x());
        }
        record_orbit(_t, h, 3600.0,satelite, obj, recoder);

        _t += h;
    }

    return;
}

fn csv_writer(data:&Vec<[f32; 10]>) -> Result<(), Box<dyn std::error::Error>> {
    let mut wtr = csv::Writer::from_path("foo.csv")?;
    
    wtr.serialize(["t", "x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2"])?;
    for iter in data {
        wtr.serialize(iter)?;
    }
    wtr.flush()?;
    Ok(())
}

fn main() {

    tester::tester();

    // let c:Vec<f64> = vec![
    //     1.0/2.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
    //     (1.0 - f64::powf(2.0, 1.0/3.0))/2.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
    //     1.0/2.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
    //     (1.0 - f64::powf(2.0, 1.0/3.0))/2.0/(2.0 - f64::powf(2.0, 1.0/3.0))
    // ];
    // let d:Vec<f64> = vec![
    //     1.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
    //     -f64::powf(2.0, 1.0/3.0)/(2.0 - f64::powf(2.0, 1.0/3.0)),
    //     1.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
    //     0.0
    // ];

    // INSTANCE_C.set(c).unwrap();
    // INSTANCE_D.set(d).unwrap();

    // // let mut pos :[f64;3] = [0.0, 0.0, 0.0];
    // // let mut _pos :[f64;3] = [0.0, 0.0, 0.0];
    // // let mut data:Vec<[f64;4]> = Vec::with_capacity(1000000);

    // let mut satelite: PhyQty::mxvr = PhyQty::mxvr {
    //     name: String::from("SATELITE"),
    //     m: 50.0,
    //     x: Vec3::Vec3([6500.0*1000.0, 0.0, 0.0]),
    //     v: Vec3::Vec3([0.0, 1000.0, 0.0]),
    //     r:0.25
    // };
    // let mut objects:Vec<PhyQty::mxvr> = Vec::new();
    // objects.push(PhyQty::mxvr { name: String::from("EARTH"), m:       0.07346*f64::powf(10.0, 24.0), x: Vec3::Vec3::new(), v: Vec3::Vec3::new(), r:  6371.0*1000.0 });
    // objects.push(PhyQty::mxvr { name: String::from("MOON") , m:       5.9724 *f64::powf(10.0, 24.0), x: Vec3::Vec3::new(), v: Vec3::Vec3::new(), r:  1737.4*1000.0 });
    // objects.push(PhyQty::mxvr { name: String::from("SUN")  , m: 1988500.0    *f64::powf(10.0, 24.0), x: Vec3::Vec3::new(), v: Vec3::Vec3::new(), r:695700.0*1000.0 });
    // let mut recorder: Vec<[f32; 10]> = Vec::with_capacity(100000);

    // cspice::call_furnsh_c();

    // sympletic4(0.1, 0.0, 100000.0, &mut satelite, &mut objects, &mut recorder);


    // if let Err(err) = csv_writer(&recorder) {
    //     println!("{}", err);
    //     std::process::exit(1);
    // }

}