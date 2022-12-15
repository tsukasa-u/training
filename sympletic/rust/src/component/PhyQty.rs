mod PhyQty {
    use crate::component::Vec3::Vec3;
    use crate::component::VecL::VecL;
    #[allow(non_camel_case_types)]
    #[allow(dead_code)]
    #[derive(Default)]
    #[derive(Clone)]
    pub struct mxvr {
        pub name: String,
        pub m: f64,
        pub x: Vec3<f64>,
        pub v: Vec3<f64>,
        pub r: f64
    }

    impl mxvr {
        #[allow(dead_code)]
        pub fn get_x(&self) -> Vec3<f64> {
            return (*self).x;
        }
        #[allow(dead_code)]
        pub fn get_v(&self) -> Vec3<f64> {
            return (*self).v;
        }
        #[allow(dead_code)]
        pub fn get_name(&self) -> &String {
            return &(*self).name;
        }
        #[allow(dead_code)]
        pub fn get_mass(&self) -> f64 {
            return (*self).m;
        }
        #[allow(dead_code)]
        pub fn get_rudius(&self) -> f64 {
            return (*self).r;
        }

        #[allow(dead_code)]
        pub fn get_x_ptr(&mut self) -> &mut Vec3<f64> {
            return &mut (*self).x;
        }
        #[allow(dead_code)]
        pub fn get_v_ptr(&mut self) -> &mut Vec3<f64> {
            return &mut (*self).v;
        }
    }

    #[allow(dead_code)]
    pub struct List<T>(Vec<T>);

    impl<T> List<T> {
        #[allow(dead_code)]
        pub fn push(&mut self, value: T) {
            (*self).0.push(value);
        }

        #[allow(dead_code)]
        pub fn Iter(&self) -> std::slice::Iter<T>{
            return (*self).0.iter();
        }
    }

    impl List<mxvr> {
        #[allow(dead_code)]
        pub fn VecL_xv_p<const N:usize>(&self) -> VecL<Vec3<f64>, N> {
            let mut tmp: VecL<Vec3<f64>, N> = VecL::new();
            for (i, v) in (*self).0.iter().enumerate() {
                // Which, Clone or Copy?
                tmp[i] = (*v).x;
                tmp[i + 1] = (*v).v;
            }
            return tmp;
        }
    }

    #[allow(unused_macros)]
    macro_rules! VecL_xp_p {
        ($x:expr) => {
            {
                let vec_size = $x.0.len();
                return $x.VecL_xv_p<vec_size>();
            }
        };
    }

}
pub use PhyQty::mxvr;