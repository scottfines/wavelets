# Discrete Wavelet Transform as a Linear Operator
===

Notes on deriving the linear operator for a given wavelt transform.

Remember from ['multi_resolution.md'] and ['cascade.md'] that we have a nice recurrence relation that is the basis of the "cascade" algorithm. 
```math
< f,\phi_{j,k} > = \sum_n h_{n-2k} < f, \phi_{j-1,k} >
```
```math
< f, \psi_{j,k} > = \sum_n (-1)^{n-2k} h_{2k-n+1} < f, \phi_{j-1,k} >
```
Suppose that we start with a set of data points $f = [f[0],f[1],...f[N-]]$ for $N=2^p$ data points. This can be written as $< f,\phi_{0,k} > = f[k]$ which 
allows us to write the recurrence relation a little more simply:
```math
f^j[k] = \sum_n h_{n-2k} f^{j-1}[k]
```
```math
d^j[k] = \sum_n (-1)^{n-2k} h_{2k-n+1} f^{j-1}[k]
```
The idea then, is to write our discrete wavelet transform as a linear matrix operation `D` where each line alternates and $f$ and a $d$ entry. Each entry
itself is determined by the properties of the wavelet in question. Because we are trying to do this on a computer, we are assuming that there are finitely 
many coefficients. Further, because our wavelet transform wants to form a basis so that we can approximate all functions, we want our individual vectors--
and therefore our matrix--should be orthonormal. That is: $\sum h_n^2 = 1$. And We redo this operation for each level, taking a smaller and smaller matrix
each time because we are removing half the entries with each iteration (the same as in the cascade approach).

We get to take advantage of another property of wavelets: the fact that $\sum_n h_h = \sqrt(2)$ (see ['matrix.md'] for why). Let us take these two equations
and do some pondering!

# Example 1: The Haar Wavelet.

Suppose that we ignore all of our other math, and just assume that we want a wavelet which has two coefficients $h_0$ and $h_1$ and is orthonormal. In that case,
we can write our recurrence relation as
```math
f^j[k] = h_0 f^{j-1}[2k] + h_1 f^{j-1}[2k+1]
```
```math
d^j[k] = h_0 f^{j-1}[2k] - h_1 f^{j-1}[2k+1]
```
(spend a bit of time with the algebra and you see how it ends up working out)
We know that $h_0^2 + h_1^2 = 1$, and that $h_0 + h_1 = \sqrt{2}$. These two equations work out to exactly one solution: $h_0 = h_1 = 1/\sqrt{2}$. This makes
our recurrence relation into
```math
f^j[k] = \frac{1}{\sqrt{2}} f^{j-1}[2k] + \frac{1}{\sqrt{2}} f^{j-1}[2k+1]
```
```math
d^j[k] = \frac{1}{\sqrt{2}} f^{j-1}[2k] - \frac{1}{\sqrt{2}} f^{j-1}[2k+1]
```
Or, written in matrix form:
```math
H = \frac{1}{\sqrt{2}}\left[ \begin{array}{cc}
1 & 1 \\
1 & -1
\end{array}
\right]
```
This matrix is written for a 2-entry $f$ array, but if you want to go to 4 entries, you would have
```math
H = \frac{1}{\sqrt{2}}\left[ \begin{array}{cc}
1 & 1 & 0 & 0 \\
1 & -1 & 0 & 0 \\
0 & 0 & 1 & 1 \\
0 & 0 & 1 & -1
\end{array}
\right]
```
And so on, where the first, 3rd,5th, and so on are the "averaging" $f^{1}$ entries and the second, 4th etc. are the $d^{1}$ entries. The discrete wavelet transform 
can thus mathematically be written as something like:

1. for $j = 0..p$ (where $N=2^p$ is the number of elements in the array) do:
  1. Form the Haar wavelet of size $2^{p-j}$
  2. Multiply $f^{j} = H * f^{j-1}$

This is mathematically nice, but not really the most effective way to do it, since matrix multiplication is typically $O(N lg N)$ and you can do the haar transform
manually in $O(N)$, but conceptually this is nice.

# Example 2: Daubechies 4-coefficient wavelet
The reason this is nice is that we can make use of this strategy to form different wavelets. Suppose that we want 4 coefficients $h_0,h_1,h_2,h_3$ in our wavelet 
decomposition. In this case, our recurrence relation gets a bit weird:
```math
f^j[k] = h_{-2k} f^{j-1}[0] + h_{1-2k} f^{j-1}[1] + h_{2-2k} f^{j-1}[2] + h_{3-2k} f^{j-1}[3]
```
```math
d^j[k] = h_{2k+1} f^{j-1}[0] +(-1)^{1-2k} h_{2k} f^{j-1}[1] + h_{2k-1} f^{j-1}[2] + (-1)^{3-2k} h_{2k-2} f^{j-1}[3]
```
Let's recall that we want a matrix operator which is orthonormal. This means that our matrix times its transpose is an identity matrix. If we look at a simple
4-element $f$ vector, we get a matrix that looks like
```math
D = \left[ \begin{array}{cc}
h_0 & h_1 & h_2 & h_3 \\
h_1 & -h_0 & h_3 & -h_2 \\
a1 & a2 & h_0 & h_1 \\
b1 & b2 & h_1 & -h_0
\end{array}
\right]
```
And we need to choose $a1,a2,b1,b2$ such that the matrix $D$ is orthogonal. This leads us to the following equations:
```math
h_0^2 + h_1^2 + h_2^2 + h_3^2 = 1
```
If we want $DD^T = I$, we thus need $h_0h_2 + h_1h_3 = 0$. So our recurrence "block" for this wavelet set is
```math
D = \left[ \begin{array}{cc}
h_0 & h_1 & h_2 & h_3 \\
h_1 & -h_0 & h_3 & -h_2 \\
h_2 & h_3 & h_0 & h_1 \\
h_3 & -h_2 & h_1 & -h_0
\end{array}
\right]
```
With the added conditions that
```math
h_0^2 + h_1^2 + h_2^2 + h_3^2 = 1
```
```math
h_0 + h_1 + h_2 + h_3 = \sqrt{2}
```
```math
h_0h_2+h_1h_3 = 0
```
Sadly, we have 3 equations and 4 unknowns, so technically we have infinitely many possibilities for the different ways that we could arrange this--that is, 
there are infinitely many wavelets that could meet this criteria (and thus be valid). 

Daubechies wavelets are specific class of this wavelet form, where an additional constraint is added: $\sum_n n*h_n = 0$. That is, that 
there are maximum "zero moments". If you work this out, you can determine the "Daubechies-4" coefficients: 
```math
h_0 = \frac{1+\sqrt{2}}{4\sqrt{2}}
```
```math
h_1 = \frac{3+\sqrt{3}}{4\sqrt{2}}
```
```math
h_2 = \frac{3-\sqrt{3}}{4\sqrt{2}}
```
```math
h_3 = \frac{1-\sqrt{3}}{4\sqrt{2}}
```

In a bit, we'll generalize this to more coefficients to form the "Daubechies Family" of wavelets, which have the 
nice property of maximum "zero moments". But first, let's look at  larger sizes of the data vector. If we imagine that $N=8$ is
the size of the vector, then we know that there are $N/2 = 4$ "avg" and "error" elements. This generates the (long) reccurrence of
```math
f^j[k] = h_{-2k} f^{j-1}[0] + h_{1-2k} f^{j-1}[1] + h_{2-2k} f^{j-1}[2] + h_{3-2k} f^{j-1}[3] + h_{4-2k} f^{j-1}[4] + h_{5-2k} f^{j-1}[5] + h_{6-2k} f^{j-1}[6] + h_{7-2k} f^{j-1}[7] 
```
```
d^j[k] = h_{2k+1} f^{j-1}[0] - h_{2k} f^{j-1}[1] + h_{2k-1} f^{j-1}[2] - h_{2k-2} f^{j-1}[3] + h_{2k-3} f^{j-1}[4] - h_{2k-4} f^{j-1}[5] + h_{2k-5} f^{j-1}[6] - h_{2k-6} f^{j-1}[7]
```
If you put in the values `0..4` into these formulas, you end up with the matrix:
```math
D = \left[ \begin{array}{cc}
h_0   & h_1   & h_2 &  h_3 &  0  &  0   & 0 & 0      \\
h_1   & -h_0  & h_3 & -h_2 &  0  &  0   & 0 & 0      \\
0     &  0    & h_0 &  h_1 & h_2 & h_3  & 0 & 0      \\
0     &  0    & h_3 & -h_2 & h_1 & -h_0 & 0 & 0      \\
0     &  0    & 0   &  0   & h_0 &  h_1 & h_2 & h_3  \\
0     &  0    & 0   &  0   & h_3 & -h_2 & h_1 & -h_0 \\
h_2   &  h_3  & 0   &  0   & 0 & 0 & h_0 &  h_1      \\
h_1   &  -h_0 & 0   &  0  & 0 & 0 & h_3 & -h_2       \\
\end{array}
\right]
```
Notice how this matrix "wraps around"--every other row moves up by two, and the last two rows pull entries from the first two values. This falls out
of the algebra, but you can also intuit that this is necessary because every value of $f$ should be used in two separate rows, but the first two
elements $f[0]$ and $f[1]$ are only used in the first rows, so you have to pick them up again at the end. This intuitition serves us well
in the actual implementation, when we can treat the last two rows as a special case and do the wrap around more efficiently.

