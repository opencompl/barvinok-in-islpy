# barvinok-in-isl
An implementation of Barvinok's Algorithm (for counting integer points in polyhedra) in islpy, preliminary to an FPL implementation.

# Documentation
## Data Structures
There are two main data structures: `Cone()` and `GenFunc()`.

`Cone` stores a signed cone using the set-of-rays representation. The `rays` attribute contains a list of the rays of the cone (all assumed to have a tail at the origin) and the `sign` attribute can take the values $\pm1$. The `d` attribute stores the dimension of the space, *i.e.*, the number of elements of each ray.  
This class has two methods:

* `get_sample_point()`: This function returns the shortest possible vector (by $l_\infty$-norm) in the space spanned by the rays of the cone. The reference for this is [the `short_vector` function in ISL](https://github.com/volcacius/barvinok/blob/cd6fd2e77b8a9e3ddc325e3636c3510c7d99bbc5/decomposer.cc#L99).
* `get_index()`: The function computes the index of the cone (to check for unimodularity). If there are more rays than dimensions of the space, it returns 0; otherwise it returns the absolute value of the determinant of the rays matrix (it is assumed that the cone is not lower-dimensional).

`GenFunc` stores a generating function as the sum of a number of terms. Mathematically, the $i^\text{th}$ term is of the form
$$\varepsilon_i \cdot \frac{\vec{x}^{\vec{v}_i}}{\prod_j (1-\vec{x}^{\vec{d}_{ij}})},$$
where $\varepsilon_i \in \{-1, 1\}$, and $\vec{v}_i, \vec{d}_{ij} \in \mathbb{R}^n$.  
The `den_exps` attribute stores the denominator exponents $\vec{d}_{ij}$, as a list of vectors for each term; the `num_exps` attribute stores the numerator exponents $\vec{v}_i$, as a vector for each term; and the `signs` attribute stores the sign of each term, as a list of $\pm1$.  
`add()` simply concatenates two `GenFunc` objects, as they both represent summations.
`__getitem__()` and `__len__()` allow us to treat a `GenFunc` object as a list of triples, consisting of a sign, a numerator exponent, and a list of denominator exponents (in that order).  

## Helper Functions
There are five helper functions, all fairly self-explanatory:

* `sign()`: Compute the sign of any nonzero integer/FP number; returns $\pm1$.
* `subAtWith()`: Substitute a list (first argument) at an index (second argument) with a value (third argument).
* `dot()`: Compute the dot product of two lists.
* `comb()`: Compute ${n \choose r}$ for any integer/FP number $n$ and any positive integer $r$.
* `conv()`: Find the convolution of two lists. If one is shorter than the other, post-pad it with zeroes and then compute the convolution.

## Main Functions
The functions are called in the following order. `count_integer_points()` takes a polytope and finds its generating function (`polytope_gen_func()`), and then finds the constant term of this function (`substitute_with_one_vector()`).

* `polytope_gen_func()`: For a given polytope, find the generating function which has one term for each point inside it. This function iterates over the vertices of the polytope (found by `cdd`'s `get_generators()` method), and at each vertex $\vec{v}$, finds the generating function, as below.
    * Find the supporting cone at $\vec{v}$ (the cone formed by the inequalities satisfied by $\vec{v}$).
    * `unimodular_decomp()`: If the index of the cone is nonzero (it is simplicial), decompose it; otherwise first triangulate it into a sum of simplicial cones, and then decompose each of them.
        * `triangulate()`: Follow the method described in Section 16.3 of [Lee & Santos, 2017](https://www.csun.edu/~ctoth/Handbook/chap16.pdf). Lifting from $\mathbb{R}^n$ to $\mathbb{R}^{n+1}$ is done by the map
        $$(v_1, \dots, v_n) \mapsto \left(v_1, \dots, v_n, \sum_{i=0}^n v_i^2\right).$$
        * `unimodular_decomp_simplicial()`: Get the shortest vector in the basis described by the cone's vectors, and replace it each ray with it by turns. If the resulting cone is full-dimensional, its index is strictly smaller than the original cone's (lower-dimensional cones can be ignored); repeat until the index is 1. This follows Algorithm 5 in [De Loera et al., 2003](https://math.ucdavis.edu/~deloera/researchsummary/barvinokalgorithm-latte1.pdf).
    * `unimodular_gen_func()`: Find the generating function of a unimodular cone (second argument) pointed at $\vec{v}$ (first argument). This is a single-term function of the type encoded in `GenFunc`, which has as exponents in the denominator the rays of the cone, and is shifted by $\vec{x}^{\vec v'}$. $\vec v'$ is defined as $\vec v$ if $\vec v$ is an integer lattice point, and by Equation (9), Section 3.1 of [Verdoolaege et al., 2007](https://link.springer.com/article/10.1007/s00453-006-1231-0) otherwise.
* Having obtained the generating function for each supporting cone, we sum them up to obtain the generating function $F$ for the polytope, as described in Theorem 2.2 of [Brion, 1988](http://www.numdam.org/article/ASENS_1988_4_21_4_653_0.pdf).
* Now, in order to find the number of terms in $F$, we substitute $\vec{x} = \vec1$. Here, we follow Section 3.2 of [Verdoolaege et al., 2007](https://link.springer.com/article/10.1007/s00453-006-1231-0).
    * `get_non_orth_vector()`: Find a vector $\vec\mu$ not orthogonal to any of the generators $\vec d_{ij}$. This is done inductively as follows:
        * If $n = 1$, return the vector $(1) \in \mathbb{R}^1$.
        * For $n > 1$, first project all the vectors $\vec d_{ij}$ down to $n-1$ dimensions $\vec d_{ij}'$. Find the vector $\vec\mu'$ not orthogonal to any $\vec d_{ij}'$. Now, we need the new coordinate $\mu_n$ such that
        $$\langle \vec\mu', \vec{d}_{ij}' \rangle + \mu_n [\vec{d}_{ij}]_n \neq 0, \forall i, j;$$
        which means that
        $$\mu_n \neq -\frac{\langle \vec\mu', \vec{d}_{ij}' \rangle}{[\vec{d}_{ij}]_n}, \forall i, j.$$
        Thus we find this set, take the maximum value in it, and increase it by one. We let $\vec\mu = (\mu'_1, \dots, \mu'_{n-1}, \mu_n)$.
    * Now, let each term in the generating function be
    $$\varepsilon \cdot \frac{\vec{x}^{\vec v}}{\prod_j (1-x^{\vec d_{j}})},$$
    where $\varepsilon \in \{-1, 1\}$, and $\vec v, \vec d_j \in \mathbb{R}^n$. We make the substitution
    $$x_i \mapsto (s+1)^{\mu_i}$$
    so that the constant term of the new expression is what we require.  

    * We have a term of the form
        $$\varepsilon * \frac{(s+1)^\text{num}}{\prod_j (1-(s+1)^{\text{den}_j})}.$$
        The numerator now has the exponent $\text{num} = \langle v, \mu \rangle$ and the denominators $\text{den}_j = \langle d_j, \mu \rangle$ for each $j$.
        * For each term $(1-(s+1)^d)$ in the denominator, if $d < 0$, we multiply the numerator and denominator by $(s+1)^{-d}$. Thus the numerator's exponent increases by $|d|$ and this term in the denominator becomes
        $$-(1-(s+1)^{|d|}),$$
        so we take the absolute values of the denominator exponents, flip the sign if necessary, and adjust the numerator exponent.
        * Now, each term in the denominator is of the form $(1-(s+1)^d)$, where $d > 0$. This is equivalent to
        $$(-s) \cdot \left(\sum_{k=0}^{d-1} (s+1)^k \right).$$
        Therefore the term now has the form
        $$\varepsilon \cdot \frac{x^\text{num}}{s^r \cdot \prod_j \left(\sum_{k=0}^{d-1} (s+1)^{k} \right)} = \varepsilon \cdot \frac{P(s)}{s^r Q(s)}.$$
        `dens` now holds $d-1$, `num` is the adjusted numerator, and `s` is the (twice) adjusted sign.
        * Finding the constant term in this expression is now reduced to finding the coefficient of $s^r$ in $\frac{P(s)}{Q(s)}$. This is done as described in Section 2.2 of [De Loera et al., 2003](https://math.ucdavis.edu/~deloera/researchsummary/barvinokalgorithm-latte1.pdf), on p. 1285.
            * We find $a_i$, the coefficients of $P$, by simply expanding the power series of $(s+1)^\text{num}$. These coefficients are then equivalent to ${\text{num} \choose i}$.
            * For the coefficients $b_i$ of $Q$, we first find the coefficients of each individual term in the product.
                * Each term in $Q$ has the form
                $$\sum_{k=0}^{\text{den}} (s+1)^k,$$
                and so the coefficient of $s^t$ in this is given by
                $$\sum_{k=0}^\text{den} {k \choose t} = {0 \choose t} + {1 \choose t} + \cdots + {\text{den} \choose t},$$
                which is equivalent to
                $${\text{den}+1 \choose k+1}.$$
                * We take the convolution over all $r$ terms to get the coefficients of $Q$.
            * Now we calculate the coefficients $c$ of $\frac{P(s)}{Q(s)}$ ($c_r$ in particular) by the recursive relation given above.
        * We sum over th constant terms of all the expressions obtained in this way.
    * This constant term is equivalent to the number of terms in the original polynomial in $x$.
* Thus we have the number of points in the polytope.