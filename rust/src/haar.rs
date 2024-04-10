// A temporary holding for the haar-based wavelet transform
// this will eventually generalize as I develop more wavelet forms for the cascades

const SQRT_2: f64 = std::f64::consts::SQRT_2;
const ROOT_2_OVER_2: f64 = 1.0 / SQRT_2;

#[derive(Debug)]
pub struct HaarWavelet {
    transform: Vec<f64>,
}

impl super::WaveletTransform for HaarWavelet {
    fn transform<T>(data: &[T]) -> Self
    where
        T: Into<f64> + Copy,
    {
        if data.len() == 0 {
            return HaarWavelet { transform: vec![] };
        }

        let floats: Vec<f64>;
        if data.len() & (data.len() - 1) != 0 {
            //the size isn't a power of two, so pad it to the nearest
            let mut p = 0;
            let mut s = 1;
            while s < data.len() {
                s <<= 1;
                p += 1;
            }

            let mut padded: Vec<f64> = Vec::with_capacity(p);
            data.into_iter().for_each(|v| padded.push((*v).into()));
            padded.resize(s, 0.0);
            floats = padded;
        } else {
            //it's already a power of 2, so we can directly transform it
            floats = data.into_iter().map(|v| (*v).into()).collect();
        }

        HaarWavelet {
            transform: dwt(&floats),
        }
    }

    fn transform_in_place(mut data: Vec<f64>) -> Self {
        if data.len() & (data.len() - 1) != 0 {
            panic!(
                "The Inverse Discrete Wavelet Transform requires that the data be a power of 2. 
                Pad out the end of the array with zero elements to ensure that this holds"
            );
        }

        dwt_in_place(&mut data);
        HaarWavelet { transform: data }
    }

    fn invert(&self) -> Vec<f64> {
        inverse_dwt(&self.transform)
    }

    fn invert_in_place(&mut self) -> Vec<f64> {
        inverse_dwt_in_place(&mut self.transform);
        std::mem::replace(&mut self.transform, vec![])
    }
}

pub fn inverse_dwt(wavelet: &[f64]) -> Vec<f64> {
    if wavelet.len() == 0 {
        return vec![]; //nothing to do
    }
    if wavelet.len() & (wavelet.len() - 1) != 0 {
        panic!(
            "The Invese Discrete Wavelet Transform requires that the data be a power of 2. 
               Pad out the end of the array with zero elements to ensure that this holds"
        );
    }
    let levels = wavelet.len().ilog2();

    let mut step: Vec<f64> = vec![];
    let mut transformed: Vec<f64> = (&wavelet[0..=1]).to_vec();
    for j in 1..=levels {
        let max_n = 1 << j;
        let last_level = max_n >> 1;
        step = Vec::with_capacity(max_n);
        for n in 0..max_n {
            let h: f64 = ROOT_2_OVER_2;
            let k: usize;
            let g: f64;
            if n % 2 == 0 {
                k = n / 2;
                g = h;
            } else {
                k = (n - 1) / 2;
                g = -h;
            }
            let c = h * transformed[k] + g * wavelet[last_level + k];
            step.push(c);
        }
        transformed = step.clone();
    }
    step
}

fn dwt(data: &[f64]) -> Vec<f64> {
    if data.len() == 0 {
        return vec![]; //nothing to do
    }
    if data.len() & (data.len() - 1) != 0 {
        panic!(
            "The Discrete Wavelet Transform requires that the data be a power of 2. 
               Pad out the end of the array with zero elements to ensure that this holds"
        );
    }
    let levels = data.len().ilog2();

    let mut diffs = Vec::with_capacity(data.len() / 2);
    diffs.resize(data.len(), 0.0);

    let mut src: Vec<f64> = data.to_vec();

    let mut split = src.len() / 2;
    for _j in 0..levels {
        for k in 0..split {
            let a: f64 = (src[2 * k] + src[2 * k + 1]) / SQRT_2;
            let d: f64 = (src[2 * k] - src[2 * k + 1]) / SQRT_2;
            src[k] = a;
            diffs[k + split] = d;
        }
        split /= 2;
    }
    diffs[0] = src[0].into();
    diffs
}

/// Compute the multi-resolution decomposition, using the Haar wavelet.
/// The Multi-resolution decomposition has lg(N)(N=2^p) entries, where each entrie
/// is a pointwise decomposition using the sums-and-diffs strategy. Each level
/// has 2^p-j average and difference terms.
fn decompose_multiresolution(data: &[f64]) -> Vec<Vec<f64>> {
    if data.len() == 0 {
        return vec![]; //nothing to do
    }
    if data.len() & (data.len() - 1) != 0 {
        panic!(
            "The Discrete Wavelet Transform requires that the data be a power of 2. 
               Pad out the end of the array with zero elements to ensure that this holds"
        );
    }

    fn level_decomp(data: &[f64]) -> Vec<f64> {
        let mut decomp = Vec::with_capacity(data.len());
        decomp.resize(data.len(), 0.0);
        let split = data.len() / 2;
        for k in 0..split {
            let a = (data[2 * k] + data[2 * k + 1]) / SQRT_2;
            let d = (data[2 * k] - data[2 * k + 1]) / SQRT_2;
            decomp[k] = a;
            decomp[split + k] = d;
        }
        decomp
    }

    let levels = data.len().ilog2() as usize;
    let mut full_decomp: Vec<Vec<f64>> = Vec::with_capacity(levels as usize);
    for j in 0..levels {
        let next_decomp: Vec<f64>;
        if j == 0 {
            next_decomp = level_decomp(data);
        } else {
            let last_decomp = &full_decomp[j - 1];
            next_decomp = level_decomp(&last_decomp[..last_decomp.len() / 2]);
        }
        full_decomp.push(next_decomp);
    }
    full_decomp
}

/// Perform the Haar Cascade wavelet transform in place.
/// This will perform the cascade algorithm to construct the Haar wavelet
/// transform on the passed in element _in place_. When the function is done,
/// the original array will be transformed into the wavelet transform. The original
/// data WILL BE LOST.
///
/// The cascade algorithm operates sequentially, where each level performs 2 steps:
/// a sums-and-diffs step to compute the values of the transform at that level, and a permutation
/// step where the sums are placed in the front, and the diffs into the back. This avoids the need
/// to have an intermediate array to hold data, but does so at the cost of (functionally) an extra
/// pass through the data, so you're trading CPU for memory. If you aren't memory constrained, it
/// is probably preferable to use a non-descructive cascade for all kinds of reasons. But
/// performance should be measured not guessed at.
fn dwt_in_place(data: &mut [f64]) {
    use crate::arrays;
    if data.len() == 0 {
        return; //nothing to do
    }
    if data.len() & (data.len() - 1) != 0 {
        panic!(
            "The Discrete Wavelet Transform requires that the data be a power of 2. 
               Pad out the end of the array with zero elements to ensure that this holds"
        );
    }

    //do a pass through the data, computing sums and differences. The sum goes into 2n, the diff
    //into 2n+1, replacing the values that were summed
    let mut len = data.len();
    while len > 1 {
        // the sum-and-diff step
        for pos in (0..len - 1).step_by(2) {
            let sum = (data[pos] + data[pos + 1]) / SQRT_2;
            let diff = (data[pos] - data[pos + 1]) / SQRT_2;
            data[pos] = sum;
            data[pos + 1] = diff;
        }
        //now the rearrangement step --only rearrange the slice that we are interested in though
        arrays::partition_evens(&mut data[0..len]);

        len /= 2;
    }
}

/// Perform the Inverse Haar Transform in place.
/// This will perform the Inverse Discrete Wavelet Transform on a haar-constructed wavelet
/// in order to reproduce the original data set. When the function is done, the
/// (wavelet) array will have been fully replaced with the inverted transform (the original data).
/// The wavelet itself WILL BE LOST.
///
/// This algorithm operates sequentially, where each level performs the inverse transform
/// of the averages and wavelet terms of the level before. A wavelet is of the form
/// [avg | wavelet terms] and is always a power of 2. Therefore, we repeatedly apply the "level
/// transform" at each level from 1 to `log(data.len())`.
///
/// The level transform is an in-place haar transform of a specific level. It operates via a
/// recursive formulation where wavelet terms are grouped together with their relevant averages,
/// and then both wavelets and averages are replaced with the inversion at that level. It works by
/// noticing that, at each level `j`, the wavelet is actually of the form
/// `[a1,a2,...,aj, w1,w2,...,wj]`, but the computation is always `ai+wi` and `ai-wi`, so we really
/// want `ai` and `wi` to be adjacent in the array. To do that, we first reorder the array so that
/// we have `[a1,a2,...aj/2,w1,w2,...,wj/2,aj/2+1,...aj,wj/2+1,...,wj]`. Then we recursively
/// reorder the left and right halves of the array until the averages are adjacent to their wavelet
/// terms, at which point we perform the calculation. In the worst case (the very last level), this
/// means doing `lg(N)` recursive reorderings, and the first reordering is of `N/2` elements, which
/// means that this algorithm runs in `Nlg(N)` time. Since the DWT also runs in Nlg(N) time, we
/// aren't _worse_ than we were before, but in principle at least a non-destructive algorithm would
/// likely prove to be faster. However, this is a case of trading time for space--if you need
/// space and have time, use this. If you need time and have space, use a non-destructive version.
fn inverse_dwt_in_place(data: &mut [f64]) {
    if data.len() == 0 {
        return; //nothing to do
    }
    if data.len() & (data.len() - 1) != 0 {
        panic!(
            "The Inverse Discrete Wavelet Transform requires that the data be a power of 2. 
               Pad out the end of the array with zero elements to ensure that this holds"
        );
    }

    fn recurse_haar(f: &mut [f64]) {
        // elements fed into this are always a power of 2 (and should always be at least a power 1
        // of 2, so there should always be at least 2 elements)
        let h = ROOT_2_OVER_2;
        if f.len() == 2 {
            let a = f[0];
            let w = f[1];
            f[0] = h * a + h * w;
            f[1] = h * a - h * w;
        } else if f.len() == 4 {
            let a1 = f[0];
            let a2 = f[1];
            let w1 = f[2];
            let w2 = f[3];
            f[0] = h * a1 + h * w1;
            f[1] = h * a1 - h * w1;
            f[2] = h * a2 + h * w2;
            f[3] = h * a2 - h * w2;
        } else {
            //we divide and conquer here. First, we split the array into quarters Q1 | Q2 | Q3 |Q4. We notice
            //that the wavelet terms for Q1 are in Q3, and the wavelet terms for Q2 are in Q4. So
            //we rearrange them to make our life easier. We rotate the middle 2 quarters so that we end up
            //with Q1 | Q3 | Q2 | Q4; This brings the wavelet terms for Q1 immediately behind the wavelet terms. From
            //there we recursively apply the haar transform on the arrays [Q1 | Q3] and [Q2 | Q4 ]
            let mid = f.len() / 2;
            let q1 = mid / 2;
            let q3 = mid + q1;
            //rotate the blocks
            f[q1..q3].rotate_right(mid - q1);

            //now recursively apply the transform on the left and right
            recurse_haar(&mut f[..mid]);
            recurse_haar(&mut f[mid..]);

            // the end result is to have an array full of averages for this level.
        }
    }

    let levels = data.len().ilog2();

    // we are counting from 1 to make the math easier to work with
    for j in 1..=levels {
        //the total number of elements at this level
        let size = 1 << j;
        recurse_haar(&mut data[..size]);
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn cascade_inverse_works() {
        let data = [1_f64, 3.0, 5.0, 11.0, 12.0, 13.0, 0.0, 1.0];

        let wavelet = super::dwt(&data);
        wavelet.iter().for_each(|k| print!("{},", k));
        println!("");

        let inverse = super::inverse_dwt(&wavelet);
        let delta = 1e-14;
        for (pos, expected) in data.iter().enumerate() {
            assert!(
                (expected - inverse[pos]).abs() < delta,
                "Element at pos {} incorrect",
                pos
            );
        }
    }

    #[test]
    fn dwt_in_place_works_simple() {
        let mut data: Vec<f64> = vec![1_f64, 3.0, 5.0, 11.0, 12.0, 13.0, 0.0, 1.0];
        let expected = data.clone();

        super::dwt_in_place(&mut data);
        data.iter().for_each(|k| print!("{},", k));
        println!("");

        let inverse = super::inverse_dwt(&data);
        let delta = 1e-14;
        for (pos, expected) in expected.iter().enumerate() {
            assert!(
                (expected - inverse[pos]).abs() < delta,
                "Element at pos {} incorrect. Expected {} but was {}",
                pos,
                expected,
                inverse[pos]
            );
        }
    }

    #[test]
    fn invert_in_place_works_simple() {
        let mut data = [1_f64, 3.0, 5.0, 11.0, 12.0, 13.0, 0.0, 1.0];
        let expected = data.clone();

        super::dwt_in_place(&mut data);
        data.iter().for_each(|k| print!("{},", k));
        println!("");

        super::inverse_dwt_in_place(&mut data);
        data.iter().for_each(|k| print!("{},", k));
        println!("");

        let delta = 1e-14;
        for (pos, expected) in expected.iter().enumerate() {
            assert!(
                (expected - data[pos]).abs() < delta,
                "Element at pos {} incorrect. Expected {} but was {}",
                pos,
                expected,
                data[pos]
            );
        }
    }

    #[test]
    fn multiresolution_decomp() {
        let data = [1_f64, 3.0, 5.0, 11.0, 12.0, 13.0, 0.0, 1.0];

        let decomp = super::decompose_multiresolution(&data);

        assert_eq!(3, decomp.len(), "Incorrect number of levels!");
        let l1: &[f64] = &[4.0, 16.0, 25.0, 1.0, -2.0, -6.0, -1.0, -1.0].map(|v| v / super::SQRT_2);
        let d1: &[f64] = &decomp[0];
        assert_eq!(l1, d1, "Incorrect decomp at level 1!");
        let l2: &[f64] = &[10.0, 13.0, -6.0, 12.0];
        let d2 = &decomp[1];
        let delta = 1e-14;
        for (pos, expected) in l2.iter().enumerate() {
            assert!(
                (expected - d2[pos]).abs() < delta,
                "Element at pos {} incorrect. Expected {} but was {}",
                pos,
                expected,
                d2[pos]
            );
        }

        let l3: &[f64] = &[23.0 / super::SQRT_2, -3.0 / super::SQRT_2];
        let d3 = &decomp[2];
        for (pos, expected) in l3.iter().enumerate() {
            assert!(
                (expected - d3[pos]).abs() < delta,
                "Element at pos {} incorrect. Expected {} but was {}",
                pos,
                expected,
                d3[pos]
            );
        }
    }
}
