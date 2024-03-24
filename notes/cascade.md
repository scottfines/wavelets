MultiResolution Analysis Algorithm
===

Here we're talking about the actual algorithm to perform a multi-resolution analysis and to get a final wavelet transform from it. In particular, we'll more-or-less manually derive the cascade algorithm.

# Recursive Computation -- The Cascade algorithm
Let's imagine that we have a function $f$, but because we're computer people, that function isn't represented continuously. Instead, we have an array of constant values $\[c_0,c_1,...\]$. For convenience,assume that there are $N = 2^p$ values in the array. We can think of this array as a sample of a continuous function $f$. Specifically, if we have a starting wavelet function $ϕ$ and corresponding wavelet basis, then we can think of these array values as $f\[n\] = \< f,\phi_{0,n} \>$ (remember that $\<\>$ is the inner product on $L^2 ( \mathbb{R} )$ )

We can think of $j$ as representing the resolution of each successive approximation, and $\< f,\phi_{j,k} \>$ as the $k$ th sample of $f$ at that resolution. That's handy! Even handier when we do a little math reduction to get

```math
< f,\phi_{j,k} > = \sum_{n} h_{n-2k} < f,\phi_{j-1,n} >
```
This is cool! It means that we can use a higher resolution approximation to build the next lower one. Recursion is the soul of algorithms, and this is perfectly recursive; we can use this to build a sequence of approximations of $f$, where each iteration has lower resolution than the previous.

Unfortunately, each of these approximations _decrease_ in resolution, which means that we are losing something in the transfer. Lets call the approximation $f^j$, and compute the difference $\delta^j = f^{j-1}-f^{j}$. It's pretty clear that $\delta^j$ captures the "error" in the $j$ th approximation step. The question just becomes how to pick out that value.

Conveniently, we happen to have this nifty function sitting around:
```math
\psi(x) = \sum_{n} (-1)^{n-1} h_{-n+1} \phi(2x-n)
```
If we get creative and define
```math
\psi_{j,k} = 2^{-j/2}\psi(2^{-j}x-k)
        = \sum g_{n-2k} \phi_{j-1,n}(x)
```
Then let's see what we get when we do the same operation here(after doing some math and using the identity of $\psi$ in terms of $\phi$):
```math
< f,\psi_{j,k} > = \sum_{n} g_{n-2k} < f,\phi_{j-1,n} >
```
If we look closely at the definitions of $g_n$ and $h_n$, we see that $\< f,𝛙_{j,k} \>$ picks out the terms that $< f,\phi_{j,k} >$ ignores. As a result, we can see that using $𝛙$ in this way allows us to pick out the error at each level of approximation.

This leads us to a fun algorithm! Suppose that we start with $f=[c_0,c_1,...]$. Assume that $< f,\phi_{0,k} > = c^0_k$ (the 0 denotes $j=0$). The algorithm is pretty straightforward:

1. compute $< f,\phi_{1,k} >$ and $< f,\psi_{1,k} >$ using the two recurrence relations above
2. repeat for $j=2,3,...$, Until we get the desired number of resolutions. 
3. Once we are done, we take our final approximation along with all of our different coefficients and that's our transform!

In principle you can stop any time you like, or could go to infinitely lower resolutions. But each time you go increase $j$, you lose half the information in the approximation. In practice, therefore, you'll eventually reach a point where there is no more information to approximate, at which point you should stop.

The basic premise here is contained in $\phi_{j,k}$. Each time we perform an iteration, we are dilating the range of $x$ that we are working over by a factor of 2. After $j=\lg(N)$ levels of approximation,our scale will be twice as large as the data set; our final approximation isn't going to change at this point because we are no longer able to take additional information into account.

Another way of looking at it is via subsampling: Each time we compute $< f,\phi_{j,k} >$, we are picking out only even terms. This means we should end up with half as many values for $f^1$ as we had for $f^0$. This will continue with each level, as each higher level of coarseness represents more data with fewer terms, until eventually there are no further terms to represent it with. This occurs at $j=lg(N)$ levels.

So the algorithm terminates after a finite number $j=N$ number of steps. Cool! But at each step, we need to compute a sum over $k$ terms. From a purely mathematical standpoint, there is no reason that $k$ cannot be infinite--we could choose a wavelet function that has infinitely many $k$ terms--but if we do, we can't compute it. Instead, we have to confine ourselves to specific wavelet functions which maintain a finite number of $h_n$ coefficients, which limits
the number of terms that are included in each subsequent operation. Thankfully, there are several wavelets that satisfy that condition. The Haar example above, for example, only has two terms for $h_n$, so the summations are quite trivial.

## Inverse Cascade Algorithm

It's really handy to be able to invert the transform, so that we can reconstruct our original data from our transform if possible. This is relatively simple when we remember that $f^{j-1} = f^j - \delta^j$. Thus, if we have each level of our coefficients, we can trivially reconstruct our original function by inverting our transforms:

```math
c^{j-1}_n =\sum_{k} h_{n-2k} c^{j}_k  + g_{n-2k} d^{j}_k
 
```
where $c^{j}_k$ is the approximation value of $f^{j}$ at location $k$, and $d^j_k$ is the error term at the same location. The algorithm is then to
compute $c^{j-1}_n$ by using the values in $c^{j}$ and $d^{j}$. After $lg(N)$ steps, you'll have $c^0_n$, which will be the original data.

That's all well and good, but that definition above is _also_ infinite. We have a finite data set, so we need to think about how we eliminate the extra terms
to build a specific inversion. 

Unfortunately, doing that requires knowing information about the specific wavelet transform that you are doing--in particular, knowing $h_n$ coefficients for that wavelet. If the wavelet is constructed with finitely many $h_n$ terms, then you can manually derive a formula for that specific wavelet transform. We can, however, build up a few intuitions which can guide our development.

If we look at the original transform, each successive approximation reduces the number of approximation _terms_ by half. So that means that the first part 
of the inverse sum $\sum h_{n-2k} c^{j}_k$ can only contain $c^{j}_k$ terms which exist at that level--half of what you're trying to build. We can thus expect that, at the $j-1$ th level of inversion, we would be generating twice as many $c^{j-1}_n$ terms as there are at the $j$ th level. This implies that our max value for $n$ should be twice what was at the previous level, which works out to be $2^{lg(N)-j}$ terms. This is nice in that it tells us how many terms we are generating for each level, but doesn't tell us anything about how many $k$ terms to sum over to generate each result, which is unfortunate. This information
is where specific knowing about the wavelet we are using becomes necessary.

### Example: Inverse Haar transform
From a previous example (#multi_resolution.md), we have a haar transform, from which we know that $h_n = 1/\sqrt{2}$ if $n = 0,1$, and $h_n = 0$ otherwise. 
Thus, we pick out only terms in the inverse in which $n-2k = 0$ and $n-2k = 1$. Recall that $n$ and $k$ are integers, which means that only _some_ of those identies are valid, depending on the value of $n$. In particular, if $n$ is even, then $k = n/2$, and if $n$ is odd, $k = (n-1)/2$ _always_. This allows us
to compute $h_{n-2k}$ and $g_{n-2k}$ trivially, depending on whether $n$ is even or odd. 

In fact, the Inverse Haar transform has a nice shortcut associated with it: It is piecewise on two elements at once. If you work out the math, you can end up with the following table for each level:
| $n$ | $k$ | $h_{n-2k}$ | $g_{n-2k}$ |
| --- | --- | --- | --- |
| 0 | 0 | $h_0$ | $h_1$ |
| 1 | 0 | $h_1$ | $-h_0$ |
| 2 | 1 | $h_0$ | $h_1$ |
| 3 | 1 | $h_1$ | $-h_0$ |

and so on, for as many terms as you need to compute your value of $n$. If $k$ is even, then $h_{n-2k} = h_0$ and $g_{n-2k} = h_1$. If $k$ is odd, then $h_{n-2k} = h_1$ and $g_{n-2k} = -h_0$. 

This makes intuitive sense. The Haar transform is taking piecewise sums and differences, so reversing that is _also_ taking piecewise sums and differences, but in a slightly different order. And since $h_0 = h_1 = 1/\sqrt{2}$, we can simplify it even further in practice.

## Example: Haar wavelets over a finite data set

Suppose $f = [1,3,5,11,12,13,0,1]$. This means that $N = 2^3=8$ values, so life is easy (we'll deal with non-dyadic arrays soon)

Using the Haar wavelet, we recall that our $h_n$ is $1/\sqrt{2}$ if $n=0,1$ and $0$ otherwise.  So $n-2k = 0,1$ or our inner product terms dissappear (thus, we only have terms for $n=2k, n=2k+1$). Similarly, our $g_n$ terms are $2k-n+1 = 0,1$
So lets' compute our transform levels:

```math
< f,ϕ_{j,k}> = \frac{1}{\sqrt{2}}(< f, ϕ_{j-1,2k} > + < f, ϕ_{j-1,2k+1} > )
```
and
```math
< f,𝛙_{j,k}> = \frac{1}{\sqrt{2}} ( < f,ϕ_{j-1,2k} > - < f,ϕ_{j-1,2k+1} > )
```
If we use our notation $c^{j}_n$ to represent the approximate value of $f^j$ at position $n$, then we can rewrite these as
```math
c^j_k = \frac{1}{\sqrt{2}}(c^{j-1}_{2k} + c^{j-1}_{2k+1} )
```
and
```math
d^j_k = \frac{1}{\sqrt{2}}(c^{j-1}_{2k} - c^{j-1}_{2k+1} )
```

So at each level, we get approximation by summing, and find out error by diffing. Sweet! Let's apply it:

```
f^0 = [1,3,5,11,12,13,0,1]
f^1 = 1/sqrt(2) [4,16,25,1]
d^1 = 1/sqrt(2) [-2,-6,-1,-1]

f^2 = 1/2 [20,26] = [10,13]
d^2 = 1/2 [-12,24] = [-6,12]

f^3 = 1/sqrt(2) [23]
d^3 = 1/sqrt(2) [-3]
```
And that's all the levels that we have. In order to represent this in a single array, we will to it tree style: first, by `j`, then increasing by `k`:
```math
haar(f) = [ 23/\sqrt{2},-3/\sqrt{2}, -6,12,-2,-6,-1,-1 ]
```

That's...handy, but how do we know that we are correct? Let's check it by computing the inverse. 

We will recall our table:

| $n$ | $k$ | $h_{n-2k}$ | $g_{n-2k}$ |
| --- | --- | --- | --- |
| 0 | 0 | $h_0$ | $h_1$ |
| 1 | 0 | $h_1$ | $-h_0$ |
| 2 | 1 | $h_0$ | $h_1$ |
| 3 | 1 | $h_1$ | $-h_0$ |

We can psuedocode up a quick algorithm for this:
```
let max_n = 2^(lg(N)-j+1)
for n in 0..<max_n:
    let k =  n even ? n/2 : (n-1)/2
    let h_{n-2k} = 1/sqrt(2)
    let g_{n-2k} = n even ? h_1 : -h_0
    c^{j-1}[n] = h_{n-2k}*c^{j}[k] + g_{n-2k}*d^{j}[k]
```

This is pretty friendly to a computer implementation, so we should be able to implement it. But let's first start doing it by hand on our previously computed wavelet transform and see if we are correct:
```math
wavelet(f) = [ 23/sqrt(2),-3/sqrt(2), -6,12,-2,-6,-1,-1 ]
```
We start with $j=2$. The maximum $n$ is thus $2$ so we have $n = 0,1$. This gives us
```math
c^2_0 = h_0 c^3_0 + h_1 d^3_0 = \frac{1}{\sqrt{2}} [ 23/\sqrt{2} -3/\sqrt{2} ] = 10
```
```math
c^2_0 = h_0 c^3_0 - h_1 d^3_0 = \frac{1}{\sqrt{2}} [ 23/\sqrt{2} +3/\sqrt{2} ] = 13
```
so $f^2 = \[ 10,13 \] $. That looks right so far! Moving to $j=1$, we get
```math
c^1_0 = 1/sqrt(2)[ c^{2}_0 + d^{2}_ ] = 1/sqrt(2) [ 10 + -6] = 4/sqrt(2)
```
```math
c^1_1 = 1/sqrt(2)[ c^{2}_0 - d^{2}_ ] = 1/sqrt(2) [ 10 + 6] = 16/sqrt(2)
```
```math
c^1_2 = 1/\sqrt{2} [ 13 + 12] = 25/sqrt(2)
```
```math
c^1_3 = 1/\sqrt{2} [ 13 - 12] = 1/sqrt(2)
```
so $f^1 = 1/\sqrt{2} \[ 4, 16, 25, 1\]$, which is also correct. The final value $j=0$ can be done as
```math
c^0_0 = 1/2[4-2] = 1
```
```math
c^0_1 = 1/2[4+2] = 3
```
```math
c^0_2 = 1/2[16 + -6] = 5
```
```math
c^0_3 = 1/2[16+6] = 11
```
```math
c^0_4 = 1/2[25-1] = 12
```
```math
c^0_5 = 1/2[25+1] = 13
```
```math
c^0_6 = 1/2[1-1] = 0
```
```math
c^0_7 = 1/2[1+1] = 1
```

so $f^0 = f = \[ 1,3,5,11,12,13,0,1 \]$ which is a perfect reconstruction! 

Thus, for the Haar wavelet at least, we have an algorithm which can perform the transform, and an approach for computing the inverse as well.


