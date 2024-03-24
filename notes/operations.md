Properties of the Wavelet Transform
===

# Storage
Given a perfect wavelet transform, the amount of storage required is $O(N)$ -- that is, a perfectly invertible wavelet transform requires as much space as the data itself. 

# Computation costs

* The cascade algorithm requires a single pass through the full data, then each subsequence pass requires half the previous data. This results in a series that looks like $N + N/2 + N/4+...2+1$ data operations. This is $O(N)$. Similarly, the reverse Cascade algorithm requires $O(N)$ operations.


# Sparsity and Compressibility.

Wavelet transforms tend to be sparse -- that is, a lot of terms are close to 0. This is pretty intuitive when you think about it--the error terms carry "difference information" where values differ from their averages. When using it on real-world data, places that are "uninteresting" are places
where two measured values are close together. Values that are close together will have small differences. So in this way, places that are "uninteresting" will have low difference terms. As an extension, dropping terms which are close to zero in the wavelet domain reduces the amount of storage
that you need without losing a lot of information. It's a form of lossy compression because _some_ information is lost, but a surprisingly low amount. This is why wavelets are a major part of modern image compression algorithms. 

So the basic data compression algorithm boils down to

0. Choose a thresholding value $\delta$
1. Perform the wavelet transform
2. Replace any difference term whose value $e \le \delta$ with 0
3. Perform the inverse transform 

The resulting data array has a lot more repeated values in it, which means that lossless compression schemes (like run-length encoding etc.) will work a lot better. There are variations on this theme of course, as people attempt to save memory during the processing and so on, but that's the basic idea.

Typically, most books and topics on Wavelets stop there--they are useful for compressing data, and that's all that we really care to talk about. But wavelet transforms are useful in other ways.

# The Error Tree
A Wavelet transform is of the form $(avg, error)$, where $avg$ is the final "average" output--the final, low-res approximation of the function, and $error$ is a list of teh difference terms that were generated during the wavelet transform. The errors are piecewise dependent on previous iterations,
so we can represent the errors as a tree.

As an intuitive example, consider $f = \[1,3,5,11,12,13,0,1\] $ and it's haar wavelet transform $haar(f) = [ 23/sqrt(2),-3/sqrt(2), -6,12,-2,-6,-1,-1 ] $. The first term in $haar(f) = 23/\sqrt{2}$ is the final
piecewise average, so write out 
```math
error(f) = [-3/\sqrt{2}, -6,12,-2,-6,-1,-1 ]
```

If you split this into the original difference levels, you get
```math
d^1 = \frac{1}{\sqrt{2}} [ -2, -6, -1, -1 ]
```
```math
d^2 = \frac{1}{2} [ -6, 12 ]
```
```math
d^3 = \frac{1}{sqrt{2}} [ -3 ]
```

If we write this as a tree (ignoring the scaling factors for now), we can see that
```
d^3[0] = -3
   |
   |- d^2[0] = -6
        |
        |- d^1[0] = -2
        |- d^1[1] = -6
   |- d^2[1] = 12
        |
        |- d^1[2] = -1
        |- d^1[3] = -1
```

If we think about values in the original space, then each node in this tree is responsible for errors that cover an interval of the original space. In face, at position `p` in the array, the corresponding
node in the tree `E[p]` will have captured difference information for the range `[p2^j,(p+1)2^j)` within the original space. 

For convenience, let's define a node in the _Error Tree_ $E$ as having the following information
```
type ErrorTreeNode:
    scale: 2^{-j/2}
    interval: [p2^j,(p+1)2^j),
    left: ErrorTree[2p+1]
    right: ErrorTree[2p+2]
```
where $j$ is determined by the number of steps through the tree that was required to reach that node.


# Answering point queries

Suppose that you have a distribution function $f$ where the positions in the array are ordered, and the values of the array are the counts for that specific position. You can use the wavelet
transform to answer the question $f[k] = ?$ for any given value of $k$ in the original space. To do this, you have the algorithm
```
func point_query(wavelet, i):
    let mut v = wavelet[0] // start with the average
    let mut j = lg(len(wavelet)) // assuming power of 2 for now because it's easier, but in practice j would start as the highest power of 2 containing the full wavelet
    let mut p = 0
    let mut sign = +1
    while p < len(wavelet) {
        let e = wavelet[1+p] // offset by 1 to account for the average term
        v += sign * (e*2^{-j/2})
        if i < (p+1)*2^{j-1} {
            //go to the left, which means the sign should be positive
            sign = +1
            p = 2*p+1 // offset by one here to account for the starting error term
        } else {
            //go to the right, which means the sign should be negative
            sign = -1
            p = 2*p+2 //offset by one to account for starting error term
        }
    }
    return v
```

This is essentially navigating the error tree, undoing the steps of the operation one by one until we reach a leaf, which we know is the original value.

This isn't _terribly_ useful, since the transform is extra work to do something that we could trivially do with the original array...unless you only have a wavelet transform. Then this is pretty nifty.
It is especially nifty if the wavelet has already been compressed--in that situation, we could compress the wavelet and _drop_ the 0 terms entirely, resulting in a smaller storage footprint while
still being able to answer point queries. 

## Answering Range queries
If you convert your original data to a cumulative distribution (such that the values contain new + all the values that came before), then you can quickly compute a range query by taking a point query
at the bottom of the range, and another at the top and subtracting the two values.

# Histograms
The property of doing point and range queries makes a wavelet very similar to a histogram, but its sparseness allows it to be more space efficient without sacrificing as much accuracy. This makes
wavelets very desirable for histogram construction.

Unfortunately, the cascade algorithm is a poor choice for histogram construction--first, it requires the same amount of storage as the data itself in order to construct (since you construct the entire
wavelet transform before truncating values). Secondly, it requires the entire data set to be organized properly into a cumulative distribution; this means an expensive pre-processing step is needed
before you can do the transform. So you can't directly perform the computation on streaming inputs. This may be tolerable for bulk processing situations (as the first step in a bulk processing pipeline maybe), but is kinda bad for real-time maintenance.

There are a few research papers that exist proposing data structures for approximate wavelet maintenance--These data structures introduce additional error in exchange for controlling the storage size
during the transform; they are probably worth exploring, but are likely not suitable for GPU computation. So we have some tradeoffs--we can have something that is suitable for processing streaming
data, or something which can be massively parallelized, but we can't have both. That's pretty par for the course with transforms really.

