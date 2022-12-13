#[macro_export]
macro_rules! concat_as_array {
    ($A:ty , $x:expr , $( $y:expr )? ; $( $z:expr ),* ) => {
        let mut ret:[A; $x + $y] = [A; Default::default()];
        let count:i32 = 0;
        if get_typename($z).starts_with("alloc::vec::Vec") {
            for iter in $z {
                ret[count] = concat_as_array!($A, $x + $y - count; $z);
                count++;
            }

        }
        if get_typename($z).starts_with("[") & get_typename($z).ends_with("[") {
            // get_typename($x).retain(|c| c != '[').retain(|c| c != ']').split(',').collect();
            for iter in $z {
                ret[count] = concat_as_array!($A, $x + $y - count; $z);
                count++;
            }
        }
        if get_typename($z) == String::from(A) {
            ret[count] = concat_as_array!($A, $x + $y - count; $z);
            count++;
        }
        return ret;
    }
}

#[allow(dead_code)]
fn get_typename<T>(_: T) -> &'static str {
    return std::any::type_name::<T>();
}