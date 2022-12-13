
pub mod mod_moons {

    // #[allow(unused_imports)]
    // use crate::component::Julian;
    // use crate::PHY_G;
    
    #[allow(non_camel_case_types)]
    pub struct orbit{
        a               :(f32, [f32; 5]),
        e               :[f32; 5],
        lambda          :(f32, [f32; 5]),
        omega_bar       :(f32, [f32; 5]),
        i               :(f32, [f32; 5]),
        capital_omega   :(f32, [f32; 5]),
        captial_gm      :f32
    }

    impl orbit {

        #[allow(dead_code)]
        pub fn new( _a               :(f32, [f32; 5]),
                    _e               :[f32; 5],
                    _lambda          :(f32, [f32; 5]),
                    _omega_bar       :(f32, [f32; 5]),
                    _i               :(f32, [f32; 5]),
                    _capital_omega   :(f32, [f32; 5]),
                    _captial_gm      :f32
        ) -> Self {
            return Self{a:_a, e:_e, lambda:_lambda, omega_bar:_omega_bar, i:_i, capital_omega:_capital_omega, captial_gm:_captial_gm};
        }

        #[allow(dead_code)]
        pub fn get_a(&self, _t:f32) -> f32 {
            return (self.a.1[0] + (self.a.1[1] + (self.a.1[2] + (self.a.1[3] + self.a.1[4]*_t)*_t)*_t)*_t)*self.a.0;
        }
        //  eccentricity
        #[allow(dead_code)]
        pub fn get_e(&self, _t:f32) -> f32 {
            return self.e[0] + (self.e[1] + (self.e[2] + (self.e[3] + self.e[4]*_t)*_t)*_t)*_t;
        }
        //  mean anomaly ? <- mean anomaly + omega + Omega where 1 >> eccentricity
        #[allow(dead_code)]
        pub fn get_lambda(&self, _t:f32) -> f32 {
            return (self.lambda.1[0] + (self.lambda.1[1] + (self.lambda.1[2] + (self.lambda.1[3] + self.lambda.1[4]*_t)*_t)*_t)*_t/3600.0)*self.lambda.0;
        }
        //  argument of pericenter ? <- omega + Omega where 1 >> inclination
        #[allow(dead_code)]
        pub fn get_omega_bar(&self, _t:f32) -> f32 {
            return (self.omega_bar.1[0] + (self.omega_bar.1[1] + (self.omega_bar.1[2] + (self.omega_bar.1[3] + self.omega_bar.1[4]*_t)*_t)*_t)*_t/3600.0)*self.omega_bar.0;
        }
        //  inclination
        #[allow(dead_code)]
        pub fn get_i(&self, _t:f32) -> f32 {
            return (self.i.1[0] + (self.i.1[1] + (self.i.1[2] + (self.i.1[3] + self.i.1[4]*_t)*_t)*_t)*_t/3600.0)*self.i.0;
        }
        //  longitude of ascending node
        #[allow(dead_code)]
        #[allow(non_snake_case)]
        pub fn get_Omega(&self, _t:f32) -> f32 {
            return (self.capital_omega.1[0] + (self.capital_omega.1[1] + (self.capital_omega.1[2] + (self.capital_omega.1[3] + self.capital_omega.1[4]*_t)*_t)*_t)*_t/3600.0)*self.capital_omega.0;
        }
        
        #[allow(dead_code)]
        #[allow(non_snake_case)]
        pub fn get_GM(&self) -> f32 {
            return self.captial_gm;
        }

        #[allow(dead_code)]
        #[allow(non_snake_case)]
        pub fn get__T(&self, _t:f32) -> f32 {
            return (self.get_a(_t).powf(3.0)*4.0*core::f32::consts::PI.powf(2.0)/self.get_GM()).powf(0.5);
        }
    }
}
pub use mod_moons::orbit;