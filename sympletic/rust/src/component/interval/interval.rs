pub mod interval {
    
    #[derive(Copy, Clone)]
    #[derive(Default)]

    
    pub struct Interval<T>(T, T);

    impl<T> Interval<T> {
        #[allow(dead_code)]
        pub fn new() -> Self where T: Default + Copy {
            return Self(Default::default(), Default::default());
        }
    }

    impl Interval<f32> {
        #[allow(dead_code)]
        pub fn reverse(self) -> Self {
            return Self(1.0/self.0, 1.0/self.1);
        }
    }
    
    impl Interval<f64> {
        #[allow(dead_code)]
        pub fn reverse(self) -> Self {
            return Self(1.0/self.0, 1.0/self.1);
        }
    }

    // impl<T> std::ops::Deref for VecN<T> {
    //     type Target = Vec<T>;
    //     fn deref(&self) -> &Self::Target {
    //         & self.0
    //     }
    // }

    // impl<T> std::ops::DerefMut for VecN<T> {
    //     fn deref_mut(&mut self) -> &mut Self::Target {
    //         &mut self.0
    //     }
    // }
    
    // impl<T> Into<VecN<T>> for Vec<T> {
    //     fn into(self) -> VecN<T> {
    //         VecN::<T>(self)
    //     }
    // }

    // // impl<T> Copy for VecN<T> where T: Copy {

    // // }

    // // impl<T> Clone for VecN<T> where T: Clone {
    // //     fn clone(&self) -> Self {
    // //         let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

    // //         let iter = self.0.iter();
    // //         for x in iter {
    // //             tmp.push((*x).clone());
    // //         }
    // //         return  tmp.into();
    // //         // return *self;
    // //     }
    // // }
    
    impl<T> std::ops::Add<Interval<T>> for Interval<T> where T: std::ops::Add<Output = T> + Copy {
        type Output = Self;
        fn add(self, rhs: Self) -> Self::Output {
            return  Interval(self.0 + rhs.0, self.1 + rhs.1);
        }
    }

    impl<T> std::ops::AddAssign<Interval<T>> for Interval<T> where T: std::ops::Add<Output = T> + Copy {
        fn add_assign(&mut self, rhs: Self) {
            *self = Interval((*self).0 + rhs.0, (*self).1 + rhs.1);
        }
    }

    impl<T> std::ops::Sub<Interval<T>> for Interval<T> where T: std::ops::Sub<Output = T> + Copy {
        type Output = Interval<T>;
        fn sub(self, rhs: Self) -> Self::Output {
            return Interval(self.0 - rhs.0, self.1 - rhs.1);
        }
    }

    impl<T> std::ops::SubAssign<Interval<T>> for Interval<T> where T: std::ops::Sub<Output = T> + Copy {
        fn sub_assign(&mut self, rhs: Self) {
            *self = Interval((*self).0 - rhs.0, (*self).1 - rhs.1);
        }
    }

    impl<T> std::ops::Neg for Interval<T> where T: std::ops::Neg<Output = T> + Copy {
        type Output = Self;

        fn neg(self) -> Self::Output {
            return  Interval(-self.0, -self.1);
        }
    }

    // impl<T> std::ops::Mul<Interval<T>> for T where T: std::ops::Mul<Output = T> {
    //     type Output = Interval<T>;
    //     fn mul(self, rhs: T) -> Self::Output {
    //         return  ;
    //     }
    // }

    impl<T> std::ops::Mul<Interval<T>> for Interval<T> where T: std::ops::Mul<Output=T> + Copy + std::cmp::PartialOrd {
        type Output = Self;
        fn mul(self, rhs: Self) -> Self::Output {
            let tmp_min: [T; 4] = [self.0*rhs.0, self.0*rhs.1, self.1*rhs.0, self.1*rhs.1];
            let tmp_max: [T; 4] = tmp_min.clone();
            return Interval(
                *tmp_min.iter().reduce(|accum, item| {
                    if accum >= item { accum } else { item }
                }).unwrap(),
                *tmp_max.iter().reduce(|accum, item| {
                    if accum <= item { accum } else { item }
                }).unwrap(),
            );
        }
    }

    impl<T> std::ops::Mul<T> for Interval<T> where T: std::ops::Mul<Output=T> + Copy {
        type Output = Interval<T>;
        fn mul(self, rhs: T) -> Self::Output {
            return  Interval(rhs*self.0, rhs*self.1);
        }
    }

    impl<T> std::ops::MulAssign<Interval<T>> for Interval<T> where T: std::ops::Mul<Output=T> + Copy + std::cmp::PartialOrd {
        fn mul_assign(&mut self, rhs: Self) {
            let tmp_min: [T; 4] = [self.0*rhs.0, self.0*rhs.1, self.1*rhs.0, self.1*rhs.1];
            let tmp_max: [T; 4] = tmp_min.clone();
            *self = Interval(
                *tmp_min.iter().reduce(|accum, item| {
                    if accum >= item { accum } else { item }
                }).unwrap(),
                *tmp_max.iter().reduce(|accum, item| {
                    if accum <= item { accum } else { item }
                }).unwrap(),
            );
        }
    }

    impl<T> std::ops::MulAssign<T> for Interval<T> where T: std::ops::Mul<Output = T> + Copy {
        fn mul_assign(&mut self, rhs: T) {
            *self = Interval(rhs*self.0, rhs*self.1);
        }
    }

    impl<T> std::ops::Div<Interval<T>> for Interval<T> where T: std::ops::Div<Output = T> + Copy + std::cmp::PartialOrd {
        type Output = Interval<T>;
        fn div(self, rhs: Self) -> Self::Output {
            let tmp_min: [T; 4] = [self.0/rhs.0, self.0/rhs.1, self.1/rhs.0, self.1/rhs.1];
            let tmp_max: [T; 4] = tmp_min.clone();
            return Interval(
                *tmp_min.iter().reduce(|accum, item| {
                    if accum >= item { accum } else { item }
                }).unwrap(),
                *tmp_max.iter().reduce(|accum, item| {
                    if accum <= item { accum } else { item }
                }).unwrap(),
            );
        }
    }

    impl<T> std::ops::Div<T> for Interval<T> where T: std::ops::Div<Output = T> + Copy {
        type Output = Interval<T>;
        fn div(self, rhs: T) -> Self::Output {
            return Self(self.0/rhs, self.1/rhs);
        }
    }
}

pub use interval::Interval;