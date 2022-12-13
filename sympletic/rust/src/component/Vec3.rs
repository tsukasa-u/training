pub mod mod_vec_3 {
    #[derive(Copy, Clone)]
    #[derive(Default)]
    pub struct Vec3<T>(pub [T; 3]) ;

    impl<T> Vec3<T> {
        pub fn as_mut_ptr(&mut self) -> &mut [T; 3] {
            return &mut (*self).0;
        }
        pub fn get(&self) -> &[T; 3] {
            return & (*self).0;
        }

        #[allow(dead_code)]
        pub fn set(&mut self, _pos: &[T; 3]) where T: Copy {
            (*self).0 = [(&_pos)[0], (&_pos)[1], (&_pos)[2]];
        }
    }

    impl Vec3<f64> {
        pub fn distance2(&self) -> f64 {
            return (*self).0[0].powf(2.0) + (*self).0[1].powf(2.0) + (*self).0[2].powf(2.0);
        }
    }

    // impl<T> Copy for Vec3<T> where T: Copy {
       
    // }

    // impl<T> Clone for Vec3<T> where T: Clone {
    //     fn clone(&self) -> Self {
    //         return *self;
    //         // Vec3([(&self).0[0].clone(), (&self).0[0].clone(), (&self).0[0].clone()])
    //     }
    // }
    

    impl<T> std::ops::Add<Vec3<T>> for Vec3<T> where T: std::ops::Add<Output = T> + Copy {
        type Output = Vec3<T>;

        fn add(self, rhs: Vec3<T>) -> Vec3<T> {
            return Self(
                [
                    self.0[0] + rhs.0[0],
                    self.0[1] + rhs.0[1],
                    self.0[2] + rhs.0[2]
                ]
            );
        }
    }

    impl<T> std::ops::AddAssign<Vec3<T>> for Vec3<T> where T: std::ops::Add<Output = T> + Copy {
        fn add_assign(&mut self, rhs: Vec3<T>) {
            *self =  Self(
                [
                    self.0[0] + rhs.0[0],
                    self.0[1] + rhs.0[1],
                    self.0[2] + rhs.0[2]
                ]
            );
        }
    }

    impl<T> std::ops::Sub<Vec3<T>> for Vec3<T> where T: std::ops::Sub<Output = T> + Copy {
        type Output = Vec3<T>;

        fn sub(self, rhs: Vec3<T>) -> Vec3<T> {
            return Self(
                [
                    self.0[0] - rhs.0[0],
                    self.0[1] - rhs.0[1],
                    self.0[2] - rhs.0[2]
                ]
            );
        }
    }

    impl<T> std::ops::SubAssign<Vec3<T>> for Vec3<T> where T: std::ops::Sub<Output = T> + Copy {
        fn sub_assign(&mut self, rhs: Vec3<T>) {
            *self =  Self(
                [
                    self.0[0] - rhs.0[0],
                    self.0[1] - rhs.0[1],
                    self.0[2] - rhs.0[2]
                ]
            );
        }
    }

    impl<T> std::ops::Neg for Vec3<T> where T: std::ops::Neg<Output = T> + Copy {
        type Output = Self;

        fn neg(self) -> Self::Output {
            return Self(
                [
                    -self.0[0],
                    -self.0[1],
                    -self.0[2]
                ]
            );
        }
    }

    // impl<T> std::ops::Mul<Vec3<T>> for T where T: std::ops::Mul<Output = T>  {
    //     type Output = Vec3<T>;

    //     fn mul(self, rhs: Vec3<T>) -> Vec3<T> {
    //         return Self(
    //             [
    //                 self * rhs.0[0],
    //                 self * rhs.0[1],
    //                 self * rhs.0[2]
    //             ]
    //         );
    //     }
    // }

    impl<T> std::ops::Mul<T> for Vec3<T> where T: std::ops::Mul<Output = T> + Copy {
        type Output = Vec3<T>;

        fn mul(self, rhs: T) -> Vec3<T> {
            return Self(
                [
                    self.0[0] * rhs,
                    self.0[1] * rhs,
                    self.0[2] * rhs
                ]
            );
        }
    }
    
    impl<T> std::ops::MulAssign<T> for Vec3<T> where T: std::ops::Mul<Output = T> + Copy {
        fn mul_assign(&mut self, rhs: T) {
            *self = Self(
                [
                    self.0[0] * rhs,
                    self.0[1] * rhs,
                    self.0[2] * rhs
                ]
            );
        }
    }

    impl<T> std::ops::Div<T> for Vec3<T> where T: std::ops::Div<Output = T> + Copy {
        type Output = Vec3<T>;

        fn div(self, rhs: T) -> Vec3<T> {
            return Self(
                [
                    self.0[0] / rhs,
                    self.0[1] / rhs,
                    self.0[2] / rhs
                ]
            );
        }
    }
}
pub use mod_vec_3::Vec3;