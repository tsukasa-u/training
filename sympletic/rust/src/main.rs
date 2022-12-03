// #[macro_use]
// extern crate lazy_static;

static INSTANCE_C: once_cell::sync::OnceCell<Vec<f32>> =  once_cell::sync::OnceCell::new();
static INSTANCE_D: once_cell::sync::OnceCell<Vec<f32>> =  once_cell::sync::OnceCell::new();

const Phy_G:f32 = 9.8;

struct Phymqp {
    m: f32,
    r: Vec3,
    v: Vec3,
    rudius: f32
}

impl Phymqp {
    fn area(&self) -> f32 {
        return self.m / self.rudius;
    }
}

#[derive(Copy, Clone)]
struct Vec3(f32, f32, f32) ;

impl std::ops::Add<Vec3> for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Vec3 {
        return Vec3(
            self.0 + rhs.0,
            self.1 + rhs.1,
            self.2 + rhs.2
        );
    }
}

impl std::ops::AddAssign<Vec3> for Vec3 {
    fn add_assign(&mut self, rhs: Vec3) {
        *self =  Self(
            self.0 + rhs.0,
            self.1 + rhs.1,
            self.2 + rhs.2
        );
    }
}

impl std::ops::Sub<Vec3> for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: Vec3) -> Vec3 {
        return Vec3(
            self.0 - rhs.0,
            self.1 - rhs.1,
            self.2 - rhs.2
        );
    }
}

impl std::ops::SubAssign<Vec3> for Vec3 {
    fn sub_assign(&mut self, rhs: Vec3) {
        *self =  Self(
            self.0 - rhs.0,
            self.1 - rhs.1,
            self.2 - rhs.2
        );
    }
}

impl std::ops::Mul<Vec3> for f32 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Vec3 {
        return Vec3(
            self * rhs.0,
            self * rhs.1,
            self * rhs.2
        );
    }
}

impl std::ops::Mul<f32> for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: f32) -> Vec3 {
        return Vec3(
            self.0 * rhs,
            self.1 * rhs,
            self.2 * rhs
        );
    }
}

fn update_orbit(_h: f32, obj:&Vec<Phymqp>) {
    for iter in obj {
        
    }
}

fn distance2(obj:&Vec3) -> f32 {
    return obj.0.powf(2.0) + obj.1.powf(2.0) + obj.2.powf(2.0);
}

fn Tp(satelite:&Phymqp) -> Vec3{
    return satelite.m*satelite.v;
}

fn Vq(satelite:&Phymqp, obj:&Vec<Phymqp>) -> Vec3{
    let mut ret:Vec3 = Vec3(0.0, 0.0, 0.0);

    for iter in obj {
        ret += (Phy_G*satelite.m/distance2(&iter.r))*iter.r;
    }

    return ret;
}

fn sympletic4(h:f32, end:f32, mut satelite:&Phymqp, obj:&Vec<Phymqp>) -> Vec<(f32, Vec<f32>, Vec<f32>)> {

    let c: &Vec<f32> = INSTANCE_C.get().unwrap();
    let d: &Vec<f32> = INSTANCE_C.get().unwrap();

    for iter in 0..4 {
        (*satelite).r += (*c)[iter]*Tp(&satelite)*h;
        (*satelite).v -= (*c)[iter]*Vq(&satelite, &obj)*h;
    }

    return vec![satelite];
}

fn main() {

    let c:Vec<f32> = vec![
        1.0/2.0/(2.0 - f32::powf(2.0, 1.0/3.0)),
        (1.0 - f32::powf(2.0, 1.0/3.0))/2.0/(2.0 - f32::powf(2.0, 1.0/3.0)),
        1.0/2.0/(2.0 - f32::powf(2.0, 1.0/3.0)),
        (1.0 - f32::powf(2.0, 1.0/3.0))/2.0/(2.0 - f32::powf(2.0, 1.0/3.0))
    ];
    let d:Vec<f32> = vec![
        1.0/(2.0 - f32::powf(2.0, 1.0/3.0)),
        -f32::powf(2.0, 1.0/3.0)/(2.0 - f32::powf(2.0, 1.0/3.0)),
        1.0/(2.0 - f32::powf(2.0, 1.0/3.0)),
        0.0
    ];

    INSTANCE_C.set(c).unwrap();
    INSTANCE_D.set(d).unwrap();

    let rect1:Phymqp = Phymqp { m: 30.0, r: Vec3(0.0, 0.0, 0.0), v:Vec3(0.0, 0.0, 0.0), rudius:10.0 };


    
    println!(
        "The area of the rectangle is {} square pixels.",
        rect1.area()
    );
}