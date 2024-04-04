Alternative approach to building Cascade algorithm
===

We start with our knowng scaling function $\phi$ from multi-resolution analysis, which gives us

```math
    \phi(x) = \sum_n h_n \phi_{j,k}(x)
```
Because we like them, we want the $\phi_{j,k}$ functions to be orthonormal, which gives us a few equations to work with:
```math
< \phi_{j,n}, \phi_{j,m} > = \int_{-\infty}&{infnty} \phi_{j,n}(x) \phi_{j,m}(x) = 1 iff m == n
```
```math
    \sum_n |h_n|^2 = 1
```
Recall that 
```math
h_n = < \phi(x), \phi_{-1,n}(x) > 
```

The first condition requires that our scaling function basis functions be orthogonal to each other, and the second requires that our coefficients are normalized so that they sum to 1. 

In general, there are allowed to be inifinitely many non-zero $h_n$ terms in this sequence. For example, the _Shannon Wavelet_ scaling function
```math
\phi(x) = \frac{sin(\pi x)}{\pi x}
```
Has inifinitely many non-zero $h_n$ coefficients. Wavelets that are of this form are sometimes referred to as _Infinite Impulse Response Filters_. 

In computational form, these sort of suck. Because you can't reduce the representation to finitely many coefficients, computing the transform of finite data involes a necessary approximation
of the function--either you do the continuous wavelet transform and use numerical integration to approximate the results, or you use the discrete wavelet transform and approximate the transform itself
with finitely many coefficients. 

But you could _also_ choose wavelets that have finitely many $h_n$ terms. It so happens that there are lots of wavelets that satisfy this--they are referred to as "Wavelets with compact support". That is,
if the Wavelet satisfies something like $\phi(x) = 0$ whenever $x$ is outside a finite range (say $[0,N]$ for some positive $N$). Some examples include
```math
\phi(x) = \begin{cases}
1 if 0 \leq x \lt 1 \\
0 otherwise
\end{cases}
```
which is the famous Haar wavelet. Another example are the Battle-Lemarie family of "spline wavelets", including
```math
\phi(x) = \begin{cases}
1-|x| if 0 \leq |x| \lt 1 \\
0 otherwise
```

Both of these trivially have compact support because they are defined to be zero outside of the interesting range. There are other famous wavelets which have this property (which we'll get to later).

So, if we are choosing a wavelet with compact support, we have finitely many $h_n$ coefficients, and we have a few useful constraints. First, we have a normality condition:
```math
\sum_n |h_n|^2 = 1
```

Secondly, we can do a little bit of fancy math on our representation of `\phi`:
```math
\int_{-\infty}^{\infty} \phi(x) dx = \sqrt{2} \sum_n h_n \int_{-\infty}^{\infty} \phi(2x-n) dx
```
if we do a neat little change of variables, and if $\int \phi \neq 0$, then we can divide off the integration and end up with
```math
\sum_n h_n = \sqrt{2}
```
Which is what I refer to as "admissibility"

If we only have 2 coefficients, then we simplify the conditions to $h_0^2 + h_1^2 = 1$ and $h_0 + h_1 = \sqrt{2}$. This is enough to uniquely define the set as $h_0 = h_1 = 1/\sqrt{2}$. Which, you'll note, is the haar wavelet. Coincidence!

But if you have 3 coefficients, then we're short an equation to solve, so we need to apply additional constraints. 

One such constraint is _invertibility_. That is, we really want the wavelet transform to be invertible. Conveniently, we enforce orthogonality, which means that
```math
    \int_{-\infty}^{\infty} \phi(x) \phi(2x-n) dx = \delta_{0,n}
```
(where $\delta_{0,m}$ is the Dirac Delta function). If you do some funky math with this, you end up with
```math
    \sum_n h_n* h_{n+2k} = \delta_{0,k}
```
(Trust me on this. The math follows pretty straightforwardly, but if you want to see the derivation, check out [1], chapter 5). Which can be referred to as the "orthogonality" condition.

Unfortunately, if you break that into 3 variables, you get
```math
    h_0^2 + h_1^2 + h_2^2 = 1
```
```math
    h_0 + h_1 + h_2 = \sqrt{2}
```
```math
    h0*h2 + h1*h3 = 0
```
That last equation spells trouble for us--$h_3 = 0$, so that means that either $h_0$ or $h_2$ is zero, and we have reduced ourselves back to the two-coefficient case. It seems that this approach will not
allow us to mechanically derive coefficient values for odd-numbers of coefficients, because that last term will always reduce the space into the next-lower even-number situation.

Doubly unfortunately, there _definitely_ exists wavelets with compact support with a non-zero integral and an odd number of coefficients: the linear Battle-Lemarie wavelet mentioned above is one such example.

So this approach has its limits. But there _are_ classes of wavelets which are very amenable to this approach. If you have an even number of coefficients, you're in good shape, although once you go to 4 coefficients, there are no more freebie equations that you can use to solve for the coefficients. for these situations, different mathematicians have adopted different approaches with different advantages and disadvantages.

Probably the most famous approach is the _Daubechies_ approach. Daubechies attempts to maximize "vanishing-moments". These are situations where a polynomial times the scaling function might have zeros within the rangeof data. It turns out that these have really nice accuracy considerations, so they are useful. The first "vanishing moment" condition is
```math
\sum_n n h_n = 0
```
You can use this to obtain the Daubechies-4(DAUB4) coefficients:
```math
h = \frac{1}{4\sqrt{2}} ( 1 + \sqrt{3}, 3+\sqrt{3}, 3-\sqrt{3}, 1-\sqrt{3} ) 
```
Applying this condition leads you to the Daubechies family of wavelets, DAUB-4,DAUB-6, etc. Using the same technique, you can derive equations for more vanishing moments (if you have $N$ coefficients, then you can have $N/2} vanishing moments) for any even-number of coefficients. The math is apparently really boring though; at this point every reference I have just points back to Daubechies' original work and just prints out a numerical table of coefficients.

So we have a technique which gives us some shortcuts to a very specific sub-families of wavelets. These wavelets are famous, but are hardly the entire possible approach. If you want a general solution, you have to go back to the multi-resolution analysis and do things the hard way. Unfortunately, for some situations, the MRA approach doesn't exactly lead to pleasant results either. The Battle-Lemarie are an example of this.

# Comments

This is a good demonstration of the difficulty in constructing a practical theory of Wavelets. You find a coherent,logical derivation (Multi-resolution analysis), only to discover that constraints that you _really_ want (like normality) aren't satisfied by mechnically applying the logic. You find a shortcut where you start from what you want--only to find that you can't get to lots of the wavelets that you want. The practical problem of how to determine the coefficients for the wavelet transform continues to vex. Every wavelet seems to have some magic math trick that makes it work out for the references, but they never show the practical end result of this. Which is probably why a lot of these wavelets are only written down mathematically and not implemented computationally anywhere. Frustrating, to say the least.

# References
[1]. Amir-Homayoon Najmi. Wavelets, a Concise Guide.
