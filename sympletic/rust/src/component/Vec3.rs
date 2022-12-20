pub mod mod_vec_3 {
    #[derive(Copy, Clone)]
    #[derive(Default)]
    #[derive(Debug)]
    pub struct Vec3<T>(pub [T; 3]) ;

    impl<T> Vec3<T> {
        pub fn new() -> Self where T: Default + Copy {
            return Self([Default::default(); 3]);
        }

        pub fn as_mut_ptr(&mut self) -> &mut [T; 3] {
            return &mut (*self).0;
        }
        #[allow(dead_code)]
        pub fn get(&self) -> &[T; 3] {
            return & (*self).0;
        }

        #[allow(dead_code)]
        pub fn set(&mut self, _pos: &[T; 3]) where T: Copy {
            (*self).0 = [(&_pos)[0], (&_pos)[1], (&_pos)[2]];
        }
    }

    impl Vec3<f64> {
        #[allow(dead_code)]
        pub fn distance(&self) -> f64 {
            return ((*self)[0].powf(2.0) + (*self)[1].powf(2.0) + (*self)[2].powf(2.0)).powf(0.5);
        }
        #[allow(dead_code)]
        pub fn distance2(&self) -> f64 {
            return (*self)[0].powf(2.0) + (*self)[1].powf(2.0) + (*self)[2].powf(2.0);
        }
        #[allow(dead_code)]
        pub fn distance3(&self) -> f64 {
            return ((*self)[0].powf(2.0) + (*self)[1].powf(2.0) + (*self)[2].powf(2.0)).powf(1.5);
        }
    }

    
    impl<T> std::ops::Deref for Vec3<T> {
        type Target = [T; 3];
        fn deref(&self) -> &Self::Target {
            return & self.0
        }
    }

    impl<T> std::ops::DerefMut for Vec3<T> {
        fn deref_mut(&mut self) -> &mut Self::Target {
            &mut self.0
        }
    }
    

    impl<T> std::ops::Add<Vec3<T>> for Vec3<T> where T: std::ops::Add<Output = T> + Copy {
        type Output = Vec3<T>;

        fn add(self, rhs: Vec3<T>) -> Vec3<T> {
            return Self(
                [
                    self[0] + rhs[0],
                    self[1] + rhs[1],
                    self[2] + rhs[2]
                ]
            );
        }
    }

    impl<T> std::ops::AddAssign<Vec3<T>> for Vec3<T> where T: std::ops::Add<Output = T> + Copy {
        fn add_assign(&mut self, rhs: Vec3<T>) {
            *self =  Self(
                [
                    self[0] + rhs[0],
                    self[1] + rhs[1],
                    self[2] + rhs[2]
                ]
            );
        }
    }

    impl<T> std::ops::Sub<Vec3<T>> for Vec3<T> where T: std::ops::Sub<Output = T> + Copy {
        type Output = Vec3<T>;

        fn sub(self, rhs: Vec3<T>) -> Vec3<T> {
            return Self(
                [
                    self[0] - rhs[0],
                    self[1] - rhs[1],
                    self[2] - rhs[2]
                ]
            );
        }
    }

    impl<T> std::ops::SubAssign<Vec3<T>> for Vec3<T> where T: std::ops::Sub<Output = T> + Copy {
        fn sub_assign(&mut self, rhs: Vec3<T>) {
            *self =  Self(
                [
                    self[0] - rhs[0],
                    self[1] - rhs[1],
                    self[2] - rhs[2]
                ]
            );
        }
    }

    impl<T> std::ops::Neg for Vec3<T> where T: std::ops::Neg<Output = T> + Copy {
        type Output = Self;

        fn neg(self) -> Self::Output {
            return Self(
                [
                    -self[0],
                    -self[1],
                    -self[2]
                ]
            );
        }
    }

    // impl<T> std::ops::Mul<Vec3<T>> for T where T: std::ops::Mul<Output = T>  {
    //     type Output = Vec3<T>;

    //     fn mul(self, rhs: Vec3<T>) -> Vec3<T> {
    //         return Self(
    //             [
    //                 self * rhs[0],
    //                 self * rhs[1],
    //                 self * rhs[2]
    //             ]
    //         );
    //     }
    // }

    impl<T> std::ops::Mul<T> for Vec3<T> where T: std::ops::Mul<Output = T> + Copy {
        type Output = Vec3<T>;

        fn mul(self, rhs: T) -> Vec3<T> {
            return Self(
                [
                    self[0] * rhs,
                    self[1] * rhs,
                    self[2] * rhs
                ]
            );
        }
    }
    
    impl<T> std::ops::MulAssign<T> for Vec3<T> where T: std::ops::Mul<Output = T> + Copy {
        fn mul_assign(&mut self, rhs: T) {
            *self = Self(
                [
                    self[0] * rhs,
                    self[1] * rhs,
                    self[2] * rhs
                ]
            );
        }
    }

    impl<T> std::ops::Div<T> for Vec3<T> where T: std::ops::Div<Output = T> + Copy {
        type Output = Vec3<T>;

        fn div(self, rhs: T) -> Vec3<T> {
            return Self(
                [
                    self[0] / rhs,
                    self[1] / rhs,
                    self[2] / rhs
                ]
            );
        }
    }
    
    impl<T> std::ops::Index<usize> for Vec3<T> {
        type Output = T;
    
        fn index(&self, _index: usize) -> &Self::Output {
            return &((*self).0[_index]);
        }
    }

    impl<T> std::ops::IndexMut<usize> for Vec3<T> {
    
        fn index_mut(&mut self, _index: usize) -> &mut Self::Output {
            return &mut ((*self).0[_index]);
        }
    }
}
pub use mod_vec_3::Vec3;