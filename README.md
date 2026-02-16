
<!-- README.md is generated from README.Rmd. Please edit that file -->

# vnorm

The goal of vnorm is to sample from variety normal distributions. For
faster computation, several models are pre-compiled upon installation.
See below for more information. The package also enables users to be
able to pre-compile models for certain polynomials and let them use it
for similar polynomials with different coefficients. Other usage include
plotting 1d varieties in 2d plots and projection onto varieties.

# Sampling from a variety

Generate samples for the variety normal distribution with mean equal to
`poly` and “standard deviation” equal to `sd`.

For a polynomial $x^2 + y^2 - 1$, we can sample using `rvnorm()`

``` r
library("vnorm")
poly <- mp("x^2 + y^2 -1")
samps <- rvnorm(2000, poly, sd = .1)

head(samps)
#>            x           y
#> 1 -1.0545771  0.05768011
#> 2  0.0340353  1.03350260
#> 3 -0.7186821 -0.69686675
#> 4  0.8587407 -0.63292490
#> 5  0.8193780 -0.60574065
#> 6  0.8092433 -0.54995186

str(samps)
#> 'data.frame':    2000 obs. of  2 variables:
#>  $ x: num  -1.055 0.034 -0.719 0.859 0.819 ...
#>  $ y: num  0.0577 1.0335 -0.6969 -0.6329 -0.6057 ...
```

Let’s plot this with ggplot:

``` r
library("ggplot2")
#> 
#> Attaching package: 'ggplot2'
#> The following object is masked from 'package:mpoly':
#> 
#>     vars
ggplot(samps, aes(x = x , y = y)) +
  geom_point() +
  coord_equal()
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

We can use pre-compiled stan models for polynomials with degree upto 3
and polynomial with up to 3 indeterminates. This helps avoid compiling
stan models.

``` r
poly <- mp("x^2 + y^2 + z^2 - 1")
samps <- rvnorm(2000, poly, sd = 0.1, pre_compiled = TRUE)
head(samps)
#>            x          y          z
#> 1  0.7337058  0.3126022  0.4641998
#> 2  1.0264048  0.4572456  0.2796481
#> 3 -0.8432177 -0.4828878  0.5819057
#> 4 -0.9009339 -0.2944216 -0.5003867
#> 5 -0.1961045 -0.7837302 -0.7070664
#> 6 -0.5652593  0.5056133 -0.7745399
```

We also want the user to be able to pre-compile models. This can be done
with `compile_stan_code`. This helps avoid re-compiling stan model for
similar polynomial with different coefficients.

``` r
poly <- mp("x^4 + y^4 - 1")
compile_stan_code(poly = poly)
#> Created registry; registered '4e9e149b8c356eb39b2171f795f8bdf4_vn.stan'
#> data {
#>   real si;
#>   real bx4;   real by4;   real b1;
#> }
#> parameters {
#>   real x;
#>   real y;
#>  }
#> model {
#>   real g = bx4*x^4+by4*y^4+b1;
#>   real dgx = 4*bx4*x^3;  real dgy = 4*by4*y^3;
#>   real ndg = sqrt(dgx^2 + dgy^2);
#>   target += normal_lpdf(0.00 | g/ndg, si); 
#> }
```

Here we have compiled stan model for polynomial of the type \$ax^4 +
by^4 -1 \$. Now, we can use the pre-compiled model with rvnorm for
similar polynomial. The `user_compiled` argument is what we will use for
this.

``` r
poly <- mp("2 x^4 + 3 y^4 - 1")
samps <- rvnorm(1000, poly = poly, sd = 0.1, user_compiled = TRUE)
head(samps)
#>            x         y
#> 1 -0.3800843 0.7903625
#> 2 -0.3893375 0.7391325
#> 3  0.8431831 0.3038976
#> 4  0.8068408 0.2984785
#> 5  0.8317430 0.3029504
#> 6  0.6251027 0.7661797
```

# Plotting

`vnorm` has `geom_variety` which is ggplot2 compliant in order to plot
1d varieties. We can plot varieties for single(mpoly object) or
multiple(mpolyList object) polynomails `geom_variety()`.

``` r
ggplot() +
 geom_variety(poly = mp("y - x^2")) +
 coord_equal()
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r
ggplot() +
  geom_variety(
    poly = mp(c("x^2 + y^2 - 1", "y - x")) ,
    xlim = c(-2, 2), ylim = c(-2, 2)
  ) +
  coord_equal()
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

## Installation

You can install the development version of vnorm from
[GitHub](https://github.com/) with:

``` r
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("dkahle/mpoly")
devtools::install_github("sonish13/vnorm)
```
