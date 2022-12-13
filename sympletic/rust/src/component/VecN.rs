pub mod mod_vec_n {
    #[derive(Clone)]
    #[derive(Default)]
    pub struct VecN<T>(Vec<T>);

    impl<T> std::ops::Deref for VecN<T> {
        type Target = Vec<T>;
        fn deref(&self) -> &Self::Target {
            & self.0
        }
    }

    impl<T> std::ops::DerefMut for VecN<T> {
        fn deref_mut(&mut self) -> &mut Self::Target {
            &mut self.0
        }
    }
    
    impl<T> Into<VecN<T>> for Vec<T> {
        fn into(self) -> VecN<T> {
            VecN::<T>(self)
        }
    }

    // impl<T> Copy for VecN<T> where T: Copy {

    // }

    // impl<T> Clone for VecN<T> where T: Clone {
    //     fn clone(&self) -> Self {
    //         let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

    //         let iter = self.0.iter();
    //         for x in iter {
    //             tmp.push((*x).clone());
    //         }
    //         return  tmp.into();
    //         // return *self;
    //     }
    // }
    
    impl<T> std::ops::Add<VecN<T>> for VecN<T> where T: std::ops::Add<Output = T> + Copy {
        type Output = VecN<T>;
        fn add(self, rhs: Self) -> Self::Output {
            let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

            let iter = self.0.iter().zip(rhs.0.iter());
            for (x, y) in iter {
                tmp.push(*x + *y);
            }
            return  tmp.into();
        }
    }

    impl<T> std::ops::AddAssign<VecN<T>> for VecN<T> where T: std::ops::Add<Output = T> + Copy {
        fn add_assign(&mut self, rhs: Self) {
            let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

            let iter = self.0.iter().zip(rhs.0.iter());
            for (x, y) in iter {
                tmp.push(*x + *y);
            }
            *self = tmp.into();
        }
    }

    impl<T> std::ops::Sub<VecN<T>> for VecN<T> where T: std::ops::Sub<Output = T> + Copy {
        type Output = VecN<T>;
        fn sub(self, rhs: Self) -> Self::Output {
            let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

            let iter = self.0.iter().zip(rhs.0.iter());
            for (x, y) in iter {
                tmp.push(*x - *y);
            }
            return  tmp.into();
        }
    }

    impl<T> std::ops::SubAssign<VecN<T>> for VecN<T> where T: std::ops::Sub<Output = T> + Copy {
        fn sub_assign(&mut self, rhs: Self) {
            let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

            let iter = self.0.iter().zip(rhs.0.iter());
            for (x, y) in iter {
                tmp.push(*x - *y);
            }
            *self = tmp.into();
        }
    }

    impl<T> std::ops::Neg for VecN<T> where T: std::ops::Neg<Output = T> + Copy {
        type Output = Self;

        fn neg(self) -> Self::Output {
            let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

            let iter = self.0.iter();
            for x in iter {
                tmp.push(-*x);
            }
            return  tmp.into();
        }
    }

    // impl<T> std::ops::Mul<VecN<T>> for T where T: std::ops::Mul<Output = T> {
    //     type Output = VecN<T>;
    //     fn mul(self, rhs: T) -> Self::Output {
    //         let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

    //         let mut iter = self.0.iter();
    //         for x in iter {
    //             tmp.push(*x * rhs);
    //         }
    //         return  tmp.into();
    //     }
    // }

    impl<T, A> std::ops::Mul<A> for VecN<T> where T: std::ops::Mul<A, Output=T> + Copy, A: Copy {
        type Output = VecN<T>;
        fn mul(self, rhs: A) -> Self::Output {
            let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

            let iter = self.0.iter();
            for x in iter {
                tmp.push(*x * rhs);
            }
            return  tmp.into();
        }
    }

    impl<T, A> std::ops::MulAssign<A> for VecN<T> where T: std::ops::Mul<A, Output = T> + Copy, A: Copy {
        fn mul_assign(&mut self, rhs: A) {
            let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

            let iter = self.0.iter();
            for x in iter {
                tmp.push(*x * rhs);
            }
            *self = tmp.into();
        }
    }

    impl<T, A> std::ops::Div<A> for VecN<T> where T: std::ops::Div<A, Output = T> + Copy, A: Copy {
        type Output = VecN<T>;
        fn div(self, rhs: A) -> Self::Output {
            let mut tmp: Vec<T> = Vec::with_capacity(self.0.len());

            let iter = self.0.iter();
            for x in iter {
                tmp.push(*x / rhs);
            }
            return  tmp.into();
        }
    }
}
pub use mod_vec_n::VecN;