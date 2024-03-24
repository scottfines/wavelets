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

