
mod component;

#[allow(unused_imports)]
use component::{Vec3, VecL, Julian};

#[allow(non_snake_case)]
mod Moon;

mod c_kernel;
use c_kernel::cspice;

// mod macros;

static INSTANCE_C: once_cell::sync::OnceCell<Vec<f64>> =  once_cell::sync::OnceCell::new();
static INSTANCE_D: once_cell::sync::OnceCell<Vec<f64>> =  once_cell::sync::OnceCell::new();

#[allow(dead_code)]
const PHY_G:f64 = 9.8;

struct Phymqp {
    name: String,
    m: f64,
    r: Vec3::Vec3<f64>,
    v: Vec3::Vec3<f64>,
    #[allow(unused)]
    rudius: f64
}

#[allow(dead_code)]
fn record_orbit(_t: f64, interval:f64, satelite: &mut Phymqp, _obj: &mut Vec<Phymqp>, recoder: &mut Vec<[f32; 10]>) {
    if _t%interval < 1.0 {
        let mut ret:[f32;10] = [0.0; 10];
        let mut iter: usize = 4;
        ret[0] = _t as f32;
        ret[1] = satelite.r.get()[0] as f32;
        ret[2] = satelite.r.get()[1] as f32;
        ret[3] = satelite.r.get()[2] as f32;
        for v in _obj {
            if v.name.contains("EARTH") == false {
                ret[iter    ] = v.r.get()[0] as f32;
                ret[iter + 1] = v.r.get()[1] as f32;
                ret[iter + 2] = v.r.get()[2] as f32;
                iter += 3;
            }
        }
        (*recoder).push(ret);
    }
}

#[allow(dead_code)]
fn update_orbit(_t: f64, _obj:&mut Vec<Phymqp>) {
    for key in _obj {
        match key.name.as_str() {
            "MOON" => {
                cspice::call_splpos_c_moon_earth_j2000((_t) as f64, key.r.as_mut_ptr());
                key.r *= 1000.0;
            },
            "SUN"  => {
                cspice::call_splpos_c_sun_earth_j2000(_t as f64, key.r.as_mut_ptr());
                key.r *= 1000.0;
            },
            _      => {}
        }
    }
}

#[allow(dead_code)]
#[allow(non_snake_case)]
fn Tp(satelite:&Phymqp) -> Vec3::Vec3<f64> {
    return satelite.v*satelite.m;
}

#[allow(dead_code)]
#[allow(non_snake_case)]
fn Vq(satelite:&Phymqp, obj:&Vec<Phymqp>) -> Vec3::Vec3<f64> {
    let mut ret:Vec3::Vec3<f64> = Vec3::Vec3([0.0, 0.0, 0.0]);

    for iter in obj {
        ret += iter.r*(PHY_G*satelite.m/iter.r.distance2());
    }

    return ret;
}

#[allow(dead_code)]
fn sympletic4(h:f64, start:f64, end:f64, satelite:&mut Phymqp, obj:&mut Vec<Phymqp>, recoder: &mut Vec<[f32; 10]>) {

    let c: &Vec<f64> = INSTANCE_C.get().unwrap();
    let d: &Vec<f64> = INSTANCE_C.get().unwrap();

    let mut _t: f64 = start;
    if end - start <= 0.0 {
        return;
    }

    while _t < end {
        for iter in 0..4 {
            (*satelite).r += Tp(&satelite)*(*c)[iter]*h;
            (*satelite).v -= Vq(&satelite, &obj)*(*d)[iter]*h;
        }

        record_orbit(_t , 3600.0,satelite, obj, recoder);

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

    let c:Vec<f64> = vec![
        1.0/2.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
        (1.0 - f64::powf(2.0, 1.0/3.0))/2.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
        1.0/2.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
        (1.0 - f64::powf(2.0, 1.0/3.0))/2.0/(2.0 - f64::powf(2.0, 1.0/3.0))
    ];
    let d:Vec<f64> = vec![
        1.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
        -f64::powf(2.0, 1.0/3.0)/(2.0 - f64::powf(2.0, 1.0/3.0)),
        1.0/(2.0 - f64::powf(2.0, 1.0/3.0)),
        0.0
    ];

    INSTANCE_C.set(c).unwrap();
    INSTANCE_D.set(d).unwrap();

    // let mut pos :[f64;3] = [0.0, 0.0, 0.0];
    // let mut _pos :[f64;3] = [0.0, 0.0, 0.0];
    // let mut data:Vec<[f64;4]> = Vec::with_capacity(1000000);


    let zeros: Vec3::Vec3<f64> = Default::default();
    let mut satelite: Phymqp = Phymqp {
        name: String::from("SATELITE"),
        m: 50.0,
        r: zeros.clone(),
        v: zeros.clone(),
        rudius:0.25
    };
    let mut objects:Vec<Phymqp> = Vec::new();
    objects.push(Phymqp { name: String::from("EARTH"), m:       0.07346*f64::powf(10.0, 24.0), r: zeros.clone(), v: zeros.clone(), rudius:  6371.0*1000.0 });
    objects.push(Phymqp { name: String::from("MOON") , m:       5.9724 *f64::powf(10.0, 24.0), r: zeros.clone(), v: zeros.clone(), rudius:  1737.4*1000.0 });
    objects.push(Phymqp { name: String::from("SUN")  , m: 1988500.0    *f64::powf(10.0, 24.0), r: zeros.clone(), v: zeros.clone(), rudius:695700.0*1000.0 });
    let mut recorder: Vec<[f32; 10]> = Vec::with_capacity(100000);

    cspice::call_furnsh_c();

    sympletic4(1.0, 0.0, 1000000.0, &mut satelite, &mut objects, &mut recorder);
    
    // for iter in 0..1000000 {
    //     cspice::call_splpos_c_moon_earth_j2000((iter) as f64, &mut pos);
    //     // cspice::call_splpos_c_sun_earth_j2000((iter*3600) as f64, &mut pos);
    //     // if (pos[0] - _pos[0]).powf(2.0) + (pos[1] - _pos[1]).powf(2.0) + (pos[2] - _pos[2]).powf(2.0) > 10000.0 {
    //     if iter%3600 == 0 {
    //         data.push([iter as f64, pos[0], pos[1], pos[2]]);
    //         _pos.clone_from(&pos);
    //     }
    // }

    if let Err(err) = csv_writer(&recorder) {
        println!("error running example: {}", err);
        std::process::exit(1);
    }



    // println!("{} {} {}", pos[0], pos[1], pos[2]);
}