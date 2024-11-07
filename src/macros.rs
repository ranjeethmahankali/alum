/// Assert that the floating point numbers are equal within the given epsilon.
#[cfg(test)]
macro_rules! assert_float_eq {
    ($a:expr, $b:expr, $eps:expr, $debug:expr) => {{
        // Make variables to avoid evaluating experssions multiple times.
        let a = $a;
        let b = $b;
        let eps = $eps;
        let error = (a - b).abs();
        if error > eps {
            eprintln!("{:?}", $debug);
        }
        assert!(
            error <= eps,
            "Assertion failed: |({}) - ({})| = {:e} <= {:e}",
            a,
            b,
            error,
            eps
        );
    }};
    ($a:expr, $b:expr, $eps:expr) => {
        $crate::macros::assert_float_eq!($a, $b, $eps, "")
    };
    ($type:ty, $a:expr, $b:expr) => {
        $crate::macros::assert_float_eq!($type, $a, $b, $type::EPSILON)
    };
}

#[cfg(test)]
macro_rules! assert_f32_eq {
    ($a:expr, $b:expr, $eps:expr, $debug:expr) => {
        $crate::macros::assert_float_eq!($a, $b, $eps, $debug)
    };
    ($a:expr, $b:expr, $eps:expr) => {
        $crate::macros::assert_float_eq!($a, $b, $eps)
    };
    ($a:expr, $b:expr) => {
        $crate::macros::assert_float_eq!($a, $b, f32::EPSILON)
    };
}

#[cfg(test)]
macro_rules! _assert_f64_eq {
    ($a:expr, $b:expr, $eps:expr, $debug:expr) => {
        $crate::macros::assert_float_eq!(f64, $a, $b, $eps, $debug)
    };
    ($a:expr, $b:expr, $eps:expr) => {
        $crate::macros::assert_float_eq!(f64, $a, $b, $eps)
    };
    ($a:expr, $b:expr) => {
        $crate::macros::assert_float_eq!(f64, $a, $b)
    };
}

#[cfg(test)]
pub(crate) use assert_f32_eq;
#[cfg(test)]
pub(crate) use assert_float_eq;
