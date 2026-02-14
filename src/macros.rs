macro_rules! nightly_exp {
    ($pub_name:ident, $inner:ident, $($args:ident),*) => {
        {
            #[cfg(feature = "nightly")]
            {
                #[cfg(feature = "std")]
                {
                    std::intrinsics::const_eval_select(($($args,)*), $inner, f64::$pub_name)
                }
                #[cfg(not(feature = "std"))]
                {
                    core::intrinsics::const_eval_select(($($args,)*), $inner, libm::$pub_name)
                }
            }
            #[cfg(not(feature = "nightly"))]
            {
                $inner($($args),*)
            }
        }
    };
}
