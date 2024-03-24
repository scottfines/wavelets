Notes on Multiresolution Analysis
===

This includes some of the basic theory of wavelets and their construction, as filtered through the lens of "how do we compute them", specifically related to Multiresolution Analysis

# An English Language description of what we are trying to do

The whole concept of "Multi-resolution analysis" is that we want to look at our data through different lenses with different resolutions. 

Imagine for a second that you're looking at a stock chart. If you are on any of those fun stock apps, they will show you a price chart with various different time-frames--last day, last week, last 30 days, last year, last 5 years, etc.
If you're looking at a chart of price over the last day, then you are doing to want some pretty good granularity in your data--something like a price taken every 5 minutes or so. But if you look over the last week, taking every 5 minutes
might be _too_ granular--you won't be able to see the "big" movements that happened over several days because you're lost in the minute-by-minute noise. You want to look at that data with a "lower resolution" --maybe every hour instead of
every 5 minutes. Similarly as you scale out to last 30 days--maybe you want only a few measurements per day now. As you scale out in time, you want to also scale out in _resolution_--lowering the resolution so that you are only picking out
the changes that take place over that time window, and ignoring the changes that take place on shorter time scales.

There are lots of algorithmic ways of accomplishing this--a popular one is to do some kind of piecewise averaging: say, taking the "highest" price over the time window and using that, or computing the average (e.g. you get an hourly reading
by taking 12 5-minute readings and averaging them together). These solutions _kinda_ work, but they tend to obscure big changes as well. If you compute an average, bigger changes within the interval will be dropped, so you may not pick up
weird flash-crashes and stuff. You've almost certainly seen this, where a given feature of a graph at one resolution just completely dissappears when you zoom out. You've probably also seen the reverse too, where trends that look really obvious
at a given scale suddenly go missing when you look at a higher resolution. What we'd really like would be a mechanism where we can zoom out to lower resolutions without losing the "big" changes that happen at the higher resolutions. Put 
another way, as we zoom out we want to keep everything _interesting_ about our previous resolution and zoom away on the uninteresting stuff so that we can really pick out what is important _across scales_.

Multi-resolution analysis is a technique for doing just that. 

# The Wavelet function

The basic idea is to break the full real space(Mathematically `L^2(ℝ)`) into "Ladder spaces" `V_i`, where 

1. $V_i \subset V_{i-1}$. That is, each successive set in the ladder gets "smaller" than the set before
2. $\bigcup V_i = L^2(ℝ)$. That is, the entire space is covered
3. $\bigcap V_i = {0}$. All sets contain 0, but that's the _only_ point covered by everything.

These are our "Resolution spaces", and they aren't really necessary except for explaining the next parts. A good example is the dyadic ranges $V_i = [-2^i,2^i]$.

So we choose a convenient "Ladder space" and then we look at specific functions $f$ which satisfy

1. $f \in V_j \iff f(2^j) \in V_0$. This more or less means that all the spaces are just scaled versions of the "central" space $V_0$. 
2. $f \in V_0 \implies f(-n) \in V_0$ for any integer $n$. That is, that linear shifts of the function stay in the space.

For example, if we choose $f_{j,k} = constant$ on the interval $[2^{j}k,2^{j}(k+1)]$,
then we can define $V_j = \\{f_{j,k} \in L^2(ℝ),k \in \mathbb{Z} \\}$, and this is a nice clean Resolution space by our definition. It's actually the "Haar space" (because of haar wavelets, but we're getting ahead of ourselves).

Finally, because we want this to be so, choose a function $\phi \in V_0$ and define
```math
\phi_{j,n}(x) = 2^{-j/2} \phi(2^{-j}x-n)
```
such that the set of $\phi_{0,n}$ forms an orthonormal basis of $V_0$. Specifically, any function $f \in V_0$ can be written as 
```math
f(x) = \sum h_n * \phi_{-1,n}(x)
```
for some $h_n ∈ \mathbb{R}$, such that $\sum |h_n|^2 = 1$. From basic linear algebra, if we have an orthonormal basis $\phi_{-1,n}$, then we can define 
```math
h_n = <\phi,\phi_{-1,n}> = \int_{-\infty}^\infty \phi(x)* \phi_{-1,n}(x) dx =\frac{1}{\sqrt{2}} \int_{-\infty}^\infty \phi(x)*\phi(2x-n)) dx
```

From here, we define a "mother wavelet function" 
```math
\psi(x) = \sum_{n\in \mathbb{Z}} (-1)^{n} h_{-n+1} * \phi(2x-n)
```
To simplify the notation, we set $g_n = (-1)^{n}h_{-n+1})$
And get
```math
\psi(x) = \sum_{n\in \mathbb{Z}} g_n \phi(2x-n)
```

So then, to find a wavelet function using the multi-resolution analysis, we start with $\phi$, then compute $h_n$; summing them up creates $\psi$. Handy!

It's worth noting that computing $h_n$ is highly dependent on the choice of $\phi$; If you choose $\phi$ poorly, there could be infinitely many $h_n$, which would lead to an incomputable wavelet transform--this might be mathematically entertaining for you, but we're aiming for some code at the end of the day here, so
we'll need to stay away from those particular wavelets.

## Example: The Haar wavelet

Start with the function $\phi(x) = 1$ if $0 \leq x < 1$, and $\phi(x) = 0$ everywhere else. This definition implies that
$\phi_{-1,n}(x)= 1$ whenever $n/2 \leq x < (n+1)/2$. So if we compute $h_n$, we get to do a bit of interval math to see that $h_n = 0$ unless $n = 0,1$, and $1/\sqrt{2}$ otherwise. Thus, we have the wavelet function

```math
\psi(x) = \frac{-1}{\sqrt{2}} \phi_{-1,0} + \frac{1}{\sqrt{2}}\phi_{-1,1}
```
working that out, we get 
```math
\psi(x) = \begin{cases}
1 & 0 \leq x \lt 1/2 \\
-1 & 1/2 \leq x \lt 1 \\
0 & otherwise
\end{cases}
```

Keep $\phi(x)$ and $\psi(x)$ in mind, because it won't be the last time we confront the haar wavelet and it's two constituent functions.

# Recursive Computation -- The Cascade algorithm
Let's imagine that we have a function `f`, but because we're computer people, that function isn't represented continuously. Instead, we have an array of constant values `[c_0,c_1,...]`. For convenience,
assume that there are `N = 2^p` values in the array. We can think of this array as a sample of a continuous function `f`. Specifically, if we have a starting wavelet function `ϕ` and corresponding
wavelet basis, then we can think of these array values as `f[n] = <f,ϕ_{0,n}>` (e.g. `j=0`). 

Now we construct the same scaling structure for `𝛙`:
```
𝛙_{j,k} = 2^{-j/2}𝛙(2^{-j}x-k)
        = 𝚺 g_{n-2k} ϕ_{j-1,n}(x)
```

Then we can think of `j` as representing the resolution of approximation, and `<f,ϕ_{j,k}>` as the `kth` sample of `f` at that resolution. That's handy! Even handier when we do a little math reduction to get
```
<f,ϕ_{j,k}> = 𝚺 h_{n-2k} < f,ϕ_{j-1,n} >
```
summing over n. This is cool! It means that we can build another approximation from the first one, of lower quality. By progressing in this way, we can build up a sequence of approximations of `f`, each with
lower and lower resolution. 

Unfortunately, each of these approximations _decrease_ in resolution, which means that we are losing something in the transfer. Lets call the approximation `f^j`, and compute the difference `𝜹^j = f^{j-1}-f^{j}`. It's pretty
clear that `𝜹^j` captures the "error" in the `j`th approximation step. The question just becomes how to pick out that value.

Conveniently, we happen to have this nifty function in mind:

```
𝛙(x) = 𝚺 (-1)^{n-1} h_{-n+1} * ϕ(2x-n)
```
We notice that our wavelet function is picking up the odd terms of the `h_n` sequence, and our approximation step is picking out the even terms. So if we define
```
𝛙_{j,k} = 2^{-j/2}𝛙(2^{-j}x-k)
```
Then let's see what we get when we do the same operation here(after doing some math and using the identity of `𝛙` in terms of `ϕ`):
```
<f,𝛙_{j,k}> = 𝚺 g_{n-2k} < f,ϕ_{j-1,n} >

```
summing over `n`. If we look closely at the definitions of `g_n` and `h_n`, we see that `<f,𝛙_{j,k}>` picks out the odd terms, while `<f,ϕ_{j,k}>` picks out the even terms. As a result, we can see that using `𝛙` in this way allows us
to pick out the error at each level of approximation.

This leads us to a fun algorithm! Suppose that we start with `f=[c_0,c_1,...]`. Assume that `<f,ϕ_{0,k}> = c_k` and call these `c^0_k` to denote being the "j=0th" coefficient. The algorithm is pretty straightforward:

1. compute `<f,ϕ_{1,k}>` and `<f,𝛙_{1,k}>` using the two recurrence relations above
2. repeat for `j=2,3,...`. Until you get the desired number of resolutions. 
3. Once we are done, we take our final approximation along with all of our different coefficients and that's our transform!

In principle you can stop any time you like. But each time you go increase `j`, you lose more information. In practice, you get to a point where the resolution is so poor that the entire function reads as zero. Once that happens,
there's no point in continuing. So where do we actually stop?

The basic premise here is contained in `ϕ_{j,k}`. Each time we perform an iteration, we are dilating the range of `x` that we are working over by a factor of 2(because we are multiplying by `2x`). After `j=N` levels of approximation,our resolution
will become twice as large as the data set, which means that our final approximation isn't going to change--from here on out, we are no longer able to add additional information to the approximation.

Another way of looking at it is via subsampling: Each time we compute `<f,ϕ_{j,k}>`, we are picking out only even terms. This means we should end up with half as many values for `f^1` as we had for `f^0`. This will continue with each level, as each
higher level of coarseness represents more data with fewer terms, until eventually there are no further terms to represent it with. This occurs at `j=N` levels.

Ok, so the algorithm terminates after a finite number `j=N` number of steps. Cool! But at each step, we need to compute a sum over `k` terms. From a purely mathematical standpoint, there is no reason that `k` cannot be infinite--we could
choose a wavelet function that has infinitely many `k` terms--but if we do, we can't compute that. Instead, we have to confine ourselves to specific wavelet functions which maintain a finite number of `h_n` coefficients, which limits
the number of terms that are included in each subsequent operation. Thankfully, there are several wavelets that satisfy that condition. The Haar example above only has two terms for `h_n`, so the summations are quite trivial.

## Inverse Cascade Algorithm

It's really handy to be able to invert the transform, so that we can reconstruct our original data from our transform if possible. This is relatively simple when we remember that `f^{j-1} = f^j - d^j`. Thus, if we have each level of our coefficients,
we can trivially reconstruct our original function by inverting our transforms:

```
<f,ϕ_{j-1,n}> =𝚺 h_{n-2k} c^{j}_k  + g_{n-2k} d^{j}_k
 
```
where $c^{j}_k$ is the approximation value of $f^{j}$ at location `k`.

## Example: Haar wavelets over a finite data set

Suppose `f = [1,3,5,11,12,13,0,1]`. This means that `N = 2^3=8` values, so life is easy (we'll deal with non-dyadic arrays soon)

Using the Haar wavelet, we recal that our `h_n` is `1/sqrt(2)` if `n=0,1` and `0` otherwise.  So `n-2k = 0,1` or our inner product terms dissappear (thus, we only have terms for `n=2k, n=2k+1`). Similarly, our `g_n` terms are `2k-n+1 = 0,1`
So lets' compute our transform levels:
```
<f,ϕ_{j,k}> = 𝚺 h_{n-2k} < f,ϕ_{j-1,n} > = h_0 <f, ϕ_{j-1,2k}> + h_1 <f, ϕ_{j-1,2k+1}>  
<f,ϕ_{j,k}> = 1/sqrt(2)*(<f, ϕ_{j-1,2k}> + <f, ϕ_{j-1,2k+1}>)

<f,𝛙_{j,k}> = 𝚺 g_{n-2k} < f,ϕ_{j-1,n} > = h_1 <f,ϕ_{j-1,2k}> - h_0<f,ϕ_{j-1,2k+1}>
<f,𝛙_{j,k}> = 1/sqrt(2) * <f,ϕ_{j-1,2k}> - <f,ϕ_{j-1,2k+1}>)
```
So at each level, we get approximation by summing, and find out error by diffing. Sweet! Let's apply it:

```
f^0 = [1,3,5,11,12,13,0,1]
f^1 =l/sqrt(2) [4,16,25,1]
d^1 = 1/sqrt(2) [-2,-6,-1,-1]

f^2 = 1/2 [20,26 ] = [10,13]
d^2 = 1/2[-12,24] = [-6,12]

f^3 = 1/sqrt(2) [ 23]
d^3 = 1/sqrt(2)[ -3]
```
And that's all the levels that we have. In order to represent this in a single array, we will to it tree style: first, by `j`, then increasing by `k`:
```
wavelet(f) = [23/sqrt(2),-3/sqrt(2), -6,12,-2,-6,-1,-1]
```

That's...handy, but how do we know that we are correct? Let's check it by computing the inverse. Let's remember that we only have `h_0(k=n/2)` and `h_1(k=n-1/2)`, giving us


```
<f,ϕ_{j-1,n}> =𝚺 h_{n-2k} c^{j}[k]  + g_{n-2k} d^{j}_k
```
where `c^{j}[k]` is the value of the approximation at level `j` at position `k`, and `d^{j}[k]` is the value of the difference term at level `j` at position `k`.

For the Haar wavelet, we have some complexity here--we have a summation over `k`, but only some values of `k` are non-zero; anywhere where `n-2k =0,1` will have a value, all others will be zero. So what values of `k` and `n` are valid for this?
We can immediately see that, for any `n`, `k = n/2` or `k=(n-1)/2`, but those two are mutually exclusive integer values--either `n/2` is an integer or `(n-1)/2` is an integer, but not both at the same time. This means that, for each `n`, we will
pick out a different `k`, creating different functions each time. 

Given that, we need to now figure out what values of `n` are valid. Remembering that we have `N` valid data points, at each level we have half the number of approximations as the level previous, and we go up `lg N` levels. This should mean that, if `k` exceeds the limit
on the number of entries in the approximation, we can terminate. That will happen after `2^{lg N-j}` entries, and `k` is determined by `n`. Thus, we should assume that `n < 2^{lgN - j}`, and that `k` is either `n/2` or `(n-1)/2` depending on whether `n` is even or odd. So
we can build out an inverse algorithm by iterating over all the valid `n` for the level, choosing a relevant `k` value, and computing the term. Written in psuedocode, we can write the algorithm for each level `j` as

```
for n in 0..<=(lg(N)-j+1):
    let k =  n even ? n/2 : (n-1)/2
    compute h_{n-2k} = 0 or 1 else continue
    compute g_{n-2k}
    c^{j-1}[n] = h_{n-2k}*c^{j}[k] + g_{n-2k}*d^{j}[k]
```

This is pretty friendly to a computer implementation, so we should be able to implement it. But let's first start doing it by hand on our previously computed wavelet transform and see if we are correct:
```
wavelet(f) = [23/sqrt(2),-3/sqrt(2), -6,12,-2,-6,-1,-1]
```

We start with `j=2`, which means `n=0,1`(`lg(8)  = 3 - 2 = 1`). For `n =0`, we have `k = 0`, giving `n-2k = 0`, and for `n=1` we have `k=0`; so we have
```
c^{2}[0] = h_0 c^{3}[0] + g_0 d^{3}[0]
         = 1/sqrt(2) [ 23/sqrt(2) + -3/sqrt(2)] = 10
c^{2}[1] = h_1 c^{3}[0] + g_1 d^{3}[0]
        = 1/sqrt(2) [ 23/sqrt(2) - -3/sqrt(2)] = 13
```
So `f^2 = [10, 13]` So far so good
For `j=1`, we get `n=0,1,2,3` (`lg(8)-1 = 2`). We have the following values for `n,k,h_{n-2k},g_{n-2k}`:
```
n = 0 => k = 0, h_0, h_1
n = 1 => k = 0, h_0, -h_0
n = 2 => k = 1, h_1, h_1
n = 3 => k = 1, h_1, -h_0
```
since `h_1` and `h_0` are the same (`1/sqrt(2)`), we can quickly write out our c values:
```
c^1[0] = 1/sqrt(2)[ c^{2}[0] + d^{2}[0] ] = 1/sqrt(2) [ 10 + -6] = 4/sqrt(2)
c^1[1] = 1/sqrt(2) [ 10 - -6] = 16/sqrt(2)
c^1[2] = 1/sqrt(2) [ 13 + 12] = 25/sqrt(2)
c^1[3] = 1/sqrt(2) [ 13 - 12] = 1/sqrt(2)

```
So `f^1 = 1/sqrt(2) [ 4,16,25,1 ]` --also lining up!

for `j=0` we have `n=0,1,2,3,4,5,6,7`, and build out the same table as
```
n = 0 => k = 0, h_0, h_1
n = 1 => k = 0, h_0, -h_0
n = 2 => k = 1, h_1, h_1
n = 3 => k = 1, h_1, -h_0
n = 4 => k = 2, h_0, h_1
n = 5 => k = 2, h_0, -h_0
n = 6 => k = 3, h_1, h_1
n = 7 => k = 3, h_1, -h_0
```
Giving
```
c^0[0] = 1/2[4-2] = 1
c^0[1] = 1/2[4+2] = 3
c^0[2] = 1/2[16 + -6] = 5
c^0[3] = 1/2[16+6] = 11
c^0[4] = 1/2[25-1] = 12
c^0[5] = 1/2[25+1] = 13
c^0[6] = 1/2[1-1] = 0
c^0[7] = 1/2[1+1] = 1
```
so `f^0 = f = [1,3,5,11,12,13,0,1]` which is a perfect reconstruction! 

Thus, for the Haar wavelet at least, we have an algorithm which can perform the transform, and an approach for computing the inverse as well.

Unfortunately, while the algorithm defines an _approach_ and a general mathematical formula, the practical implementation is unique to each choice of wavelet because you need to exploit the mathematical properties of each
wavelet that you select--the coefficients of the transform `h_n` have to be computed ahead of time, and the inverse algorithm in particular is sensitive to the wavelet's nature. So math doesn't save us entirely, although
it does create a space for us.

### 1-Dimensional Haar Wavelet as a matrix transform

The Haar wavelet pointwise sums two terms, and then pointwise differences them. For a vector of two variables `f = [a[0],a[1]]`, this looks like the matrix 
```
H = 
    | 1  1 |
    | 1 -1 |

```

You can manually verify that this does the piecewise summation and differencing correctly (and we scale the final result by `1/sqrt(2)` to maintain orthonormality).

If `f` is a vector of `4` elements, then we can repeat this by having a block-wise matrix that looks like

```
H = 
    | 1  1  0  0 |
    | 1 -1  0  0 |
    | 0  0  1  1 |
    | 0  0  1 -1 |

```
Which will produce the pointwise sums and differences, but not in a particularly friendly order--the approximations and differences are mixed up. Rearranging the rows a bit to get

```
H = 
    | 1  1  0  0 |
    | 0  0  1  1 |
    | 1 -1  0  0 |
    | 0  0  1 -1 |

```
In general, we can scale this to any size `N=2^p` by doing the piecewise formula
```
H = 
    | 1  1  0  0  ...|
    | 0  0  1  1  ...|
    | 0  0  0  0  1 1 ...|
    | ... |
    | 1 -1  0  0  ...|
    | 0  0  1 -1  ...|
    | 0  0  0  0  1 -1 ...|

```

Will group the approximations and error terms together in the same way that the cascade algorithm will present them. However, at each level, the number of elements in the resulting vector which are used in the transform decrease by half, so each iteration of the transform
will use the Haar matrix that is half the size. After `j` iterations, there will be a `2^{-j/2}` scaling factor in the front and the haar transform as a vector.

This approach has the extremely pleasant property of just being a sequence of matrix multiplications, which has _many_ excellent algorithms and implementations for extra efficiency. 


So we have two algorithms--one is the perform the cascade approach and scale up the approximations at each level. The other is to form the relevant Wavelet transform matrix, and repeatedly apply it to values as a vector. That's pretty awesome, although the matrix multiplication
approach requires us to hand-build the correct matrix of the correct size each time, which may require more memory than we'd want to do, and in practice the cascades algorithm should have the same piecewise-parallelism capabilities as the matrix multiplication, so we should be able
to create a parallel cascade algorithm in a relatively straightforward manner. We'll have to do that after we implement the raw cascade multi-level resolution approach.

# Commentary
I think MultiResolution Analysis is actually a lot more interesting than people give it credit. Most of the mathematical texts only talk about it as a means to an end--use MRA to get to the wavelet transform, and then talk about the properties of the transform itself. But the idea
that we can zoom in and out dynamically across the same data without losing interesting elements is a compelling idea to me. Imagine being able to essentially layer on the resolutions, starting with a lower resolution and increasing it _only_ in the parts that are interesting to you,
and being able to go back and forth without needing to recompute tons of data points? That's a cool idea, and I can see how it would make for neat visualizations of lots of data--particularly sensor data. You can start to poke out trends _and_ deviations from that trend at the same time. 
It's probably limited as an _analytical_ tool--I doubt that using the intermediate steps of the MRA would allow you to speed up any particular computations when doing (say) compression or machine learning or whatever, but it would have a powerful _explanatory_ effect, which is something
that I definitely think is missing from a lot of statistical techniques that are currently floating around.

Another problem with MRA that I find is that there really isn't a good explanation of what it's doing--it's either a 1 paragraph note for journalists which explains nothing, or it's a phD thesis complete with proofs and Hilbert spaces. There's no really clear, intuitive explanation for
the process. And of course, since there are other algorithms for computing wavelets (lifting schemes) and the choice of wavelet matters, they tend to fall into the rabbit hole of over-generalization or over-specialization.

# References
[1] Daubechies. _10 Lectures on Wavelets_.

