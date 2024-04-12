// internal modules
mod arrays;
mod haar;

/// A Wavelet Transform.
///
pub trait WaveletTransform {
    /// Perform the Discreate Wavelet Transform on the specified data non-destructively.
    ///
    /// This implementation should create and manage its own memory, leaving 'data' unaffected.
    /// However, be careful of the larger memory footprint for larger data set in particular.
    fn transform<T>(data: &[T]) -> Self
    where
        T: Into<f64> + Copy;

    /// Perform the Wavelet Transform on the passed-in data in place destructively.
    ///
    /// The returned type should retain ownership of the 'data' object throughout its
    /// lifetime.
    fn transform_in_place(data: Vec<f64>) -> Self;

    /// Invert the Transform.
    ///
    /// This operation is performed on a copy of the data, which does not destroy this instance.
    ///
    /// The returned vector will hold the original data as reconstructed by this wavelet
    /// (modulo errors introduced by floating point operations).
    fn invert(&self) -> Vec<f64>;

    /// Invert the Transform in place.
    ///
    /// This operations is performed on the underlying data in place, and the underlying memory
    /// is moved directly into the return value (destroying the viability of the transform object
    /// itself).
    ///
    /// Note that in-place inversions are typically more computationally expensive than copy-based
    /// approaches, since they restrict the total amount of memory used during the operation.
    /// So use this if you are memory-constrained but willing to deal with additional CPU steps. If
    /// you are not in that situation, you are likely better off using [`invert`] instead.
    fn invert_in_place(&mut self) -> Vec<f64>;
}

/// A Multi-Resolution Decomposition.
///
/// A Multi-Resolution Decomposition has all the levels of the discrete wavelet transform
/// available for analysis, and can be converted into a type of WaveletTransform.
pub trait MRDecomposition<WT: WaveletTransform>: Into<WT> {
    fn decompose<T>(data: &[T]) -> Self
    where
        T: Into<f64> + Copy;
}

/// Perform the Discrete Wavelet Transform(DWT) on the specified data.
///
/// This function performs the transform in a non-destructive way, creating a new
/// data structure to hold the underlying data. This function is preferable when memory is not
/// a significant concern, as it tends to be algorithmically more efficient than an in-place
/// solution. However, if you are afraid about memory usage (for example, if you are transforming
/// very large data sets or operating on a device with very limited memory capabilities), you may
/// prefer [`dwt_in_place`], as it will reuse an existing memory space to perform its underlying
/// calculation.
///
/// Note that the discrete wavelet transform requires that the input data be a power of 2. This
/// function deals with this by padding--if the data isn't sized to be a power of two, then a new
/// vector will be created which is the next power of 2 higher, and the result will be padded with
/// zeros. This padding doesn't significantly affect the resulting transform, but ensures that the
/// calculations can be correctly executed.
pub fn dwt<T, W>(data: &[T]) -> W
where
    T: Into<f64> + Copy,
    W: WaveletTransform,
{
    WaveletTransform::transform(data)
}

/// Perform the Discrete Wavelet Transform(DWT) on the specified data in place.
///
/// This function produces an identical result as [`dwt`], but performs all of its operations
/// in place using the memory which is passed (taking ownership of that memory in the process).
/// This is suitable when memory is a more significant resource than CPU, because most in place
/// transform implementations must perform additional steps to correctly perform the
/// transform/inverse, which results in poorer algorithmic performance than using a copy of the
/// memory space.
///
/// Note that the discrete wavelet transform requires that the input data be a power of 2. Because
/// this function operate in place without allocating new memory, this function requires that the
/// data already be a power of 2--otherwise, the function will panic.
pub fn dwt_in_place<W>(mut data: Vec<f64>) -> W
where
    W: WaveletTransform,
{
    if data.len() & (data.len() - 1) != 0 {
        panic!(
            "The Inverse Discrete Wavelet Transform requires that the data be a power of 2. 
               Pad out the end of the array with zero elements to ensure that this holds"
        );
    }

    WaveletTransform::transform_in_place(data)
}
