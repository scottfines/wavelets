/// Utilities for manipulating arrays, used in the wavelet transforms.

/// Permute the elements in the slice such that even-numbered elements
/// are moved to the front, and odd-numbered elements are moved to the back,
/// while otherwise maintaining their original sort order.
/// For example, if you start with [0,1,2,3,4,5], this will sort the array
/// into [0,2,4,1,3,5].
/// This tranformation is performed _in place_, so the input slice will not
/// have elements in the same location when this function call is done.
pub fn partition_evens<E>(d: &mut [E])
where
    E: Default,
{
    use std::mem;
    if d.len() < 2 {
        //nothing to do
        return;
    }

    let mut p = 1;
    let mut n = 0;
    let mid = if d.len() % 2 == 0 {
        d.len() / 2
    } else {
        d.len() / 2 + 1
    };

    while p < mid {
        let sn = n;
        let mut t = mem::replace(&mut d[p], E::default());
        loop {
            //compute destination
            let dp = if p % 2 == 0 { n } else { mid + n };
            let t2 = mem::replace(&mut d[dp], E::default());
            d[dp] = t;
            t = t2;
            p = dp;
            n = if dp % 2 == 0 { dp / 2 } else { (dp - 1) / 2 };

            if n == sn {
                break;
            }
        }
        n += 1;
        p = 2 * n + 1;
    }
}

pub fn pad<T>(d: &[T]) -> Vec<f64>
where
    T: Into<f64> + Copy,
{
    if d.len() & (d.len() - 1) != 0 {
        //the size isn't a power of two, so pad it to the nearest
        let mut p = 0;
        let mut s = 1;
        while s < d.len() {
            s <<= 1;
            p += 1;
        }

        let mut padded: Vec<f64> = Vec::with_capacity(p);
        d.into_iter().for_each(|v| padded.push((*v).into()));
        padded.resize(s, 0.0);
        padded
    } else {
        d.into_iter().map(|v| (*v).into()).collect()
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn partition_works_even() {
        let mut data = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];

        super::partition_evens(&mut data);

        assert_eq!(
            [0.0, 2.0, 4.0, 6.0, 1.0, 3.0, 5.0, 7.0],
            data,
            "Incorrect permutation!"
        );
    }

    #[test]
    fn partition_works_odd() {
        let mut data = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];

        super::partition_evens(&mut data);

        assert_eq!(
            [0.0, 2.0, 4.0, 6.0, 8.0, 1.0, 3.0, 5.0, 7.0],
            data,
            "Incorrect permutation!"
        );
    }

    #[test]
    fn partition_size_3() {
        let mut data = [0.0, 1.0, 2.0];

        super::partition_evens(&mut data);

        assert_eq!([0.0, 2.0, 1.0], data, "Incorrect permutation!");
    }

    #[test]
    fn partition_several() {
        for size in 0..10 {
            let mut data: Vec<usize> = Vec::new();
            for i in 0..=size {
                data.push(i);
            }

            super::partition_evens(&mut data[..]);

            //the first size/2 elements should be even then odd
            for i in 0..=size / 2 {
                assert!(
                    data[i] % 2 == 0,
                    "Element d[{}] should be even, but is {}!",
                    i,
                    data[i]
                );
            }
            for i in size / 2 + 1..size {
                assert!(
                    data[i] % 2 != 0,
                    "Element d[{}] should be odd, but is {}!",
                    i,
                    data[i]
                );
            }
        }
    }
}
