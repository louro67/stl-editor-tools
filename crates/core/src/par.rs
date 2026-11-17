//! Data parallelism when the `parallel` feature is on, plain loops otherwise.

#[cfg(feature = "parallel")]
pub fn map<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
    use rayon::prelude::*;
    (0..n).into_par_iter().map(f).collect()
}

#[cfg(not(feature = "parallel"))]
pub fn map<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync + Send) -> Vec<T> {
    (0..n).map(f).collect()
}

/// Apply `f` to each chunk of `v` (chunk index, chunk).
#[cfg(feature = "parallel")]
pub fn chunks_mut<T: Send>(v: &mut [T], size: usize, f: impl Fn(usize, &mut [T]) + Sync + Send) {
    use rayon::prelude::*;
    v.par_chunks_mut(size)
        .enumerate()
        .for_each(|(i, c)| f(i, c));
}

#[cfg(not(feature = "parallel"))]
pub fn chunks_mut<T: Send>(v: &mut [T], size: usize, f: impl Fn(usize, &mut [T]) + Sync + Send) {
    v.chunks_mut(size).enumerate().for_each(|(i, c)| f(i, c));
}
