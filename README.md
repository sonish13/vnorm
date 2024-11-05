
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

For a polynomial $x^2 + y^2 -1$, we can sample using `rvnorm`

``` r
library(vnorm)
poly <- mpoly::mp("x^2 + y^2 -1")
samps <- rvnorm(2000, poly, sd = .1)
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 finished in 0.2 seconds.
#> Chain 2 finished in 0.2 seconds.
#> Chain 3 finished in 0.3 seconds.
#> Chain 4 finished in 0.2 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.2 seconds.
#> Total execution time: 1.2 seconds.
head(samps)
#>           x         y
#> 1  0.966994  0.405642
#> 2  0.958511  0.375407
#> 3  0.394862 -0.967098
#> 4 -0.736241  0.734586
#> 5 -0.732316  0.734969
#> 6 -0.707225  0.701170
str(samps)
#> 'data.frame':    2000 obs. of  2 variables:
#>  $ x: num  0.967 0.959 0.395 -0.736 -0.732 ...
#>  $ y: num  0.406 0.375 -0.967 0.735 0.735 ...
```

Let’s plot this with ggplot:

``` r
library(ggplot2)
#> 
#> Attaching package: 'ggplot2'
#> The following object is masked from 'package:mpoly':
#> 
#>     vars
samps |> 
  ggplot(aes(x = x , y = y)) +
  geom_point() +
  coord_equal()
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

We can use pre-compiled stan models for polynomials with degree upto 3
and polynomial with upto 3 indeterminates. This helps avoid compiling
stan models.

``` r
poly <- mpoly::mp("x^2 + y^2 + z^2 - 1")
samps <- rvnorm(2000, poly, sd = 0.1, pre_compiled = TRUE)
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 finished in 0.3 seconds.
#> Chain 2 finished in 0.4 seconds.
#> Chain 3 finished in 0.3 seconds.
#> Chain 4 finished in 0.3 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.3 seconds.
#> Total execution time: 1.6 seconds.
head(samps)
#>           x         y          z
#> 1  0.162459 -0.988889 -0.0092164
#> 2 -0.497968 -0.547431  0.6876140
#> 3  0.352122 -0.782564  0.5125490
#> 4 -0.470062  0.781425  0.4553070
#> 5 -0.055613  0.521399 -0.8016060
#> 6  0.126799  0.060375 -0.9487200
```

We also want the user to be able to pre-compile models. This can be done
with `compile_stan_code`. This helps avoid re-compiling stan model for
similar polynomial with different coefficients.

``` r
poly <- mpoly::mp("x^4 + y^4 - 1")
compile_stan_code(poly = poly)
#> compiled_stan_info variable created in global environment
#> [1] "Model Compiled"
```

Here we have compiled stan model for polynomial of the type \$ax^4 +
by^4 -1 \$. Now, we can use the pre-compiled model with rvnorm for
similar polynomial. The `user_compiled` argument is what we will use for
this.

``` r
poly <- mpoly::mp("2 x^4 + 3 y^4 - 1")
samps <- rvnorm(1000, poly = poly, sd = 0.1, user_compiled = TRUE)
#> Running MCMC with 4 sequential chains...
#> 
#> Chain 1 finished in 0.1 seconds.
#> Chain 2 finished in 0.1 seconds.
#> Chain 3 finished in 0.1 seconds.
#> Chain 4 finished in 0.0 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.1 seconds.
#> Total execution time: 0.6 seconds.
head(samps)
#>          x          y
#> 1 0.861923  0.3347540
#> 2 1.010160  0.3768460
#> 3 0.744661  0.0496880
#> 4 0.993849  0.0955615
#> 5 0.961165  0.0962827
#> 6 0.937857 -0.2780060
```

# Plotting

`vnorm` has `geom_variety` which is ggplot2 compliant in order to plot
1d varieties. We can plot varieties for single(mpoly object) or
multiple(mpolyList object) polynomails `geom_variety()`.

``` r
poly <- mpoly::mp("y - x^2")
ggplot() +
 geom_variety(poly = poly) +
 coord_equal()
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r
p1 <- mp("x^2 + y^2 - 1")
p2 <- mp("y - x")
poly <- mpolyList(p1, p2)
ggplot() +
  geom_variety(poly = poly , xlim = c(-2, 2), ylim = c(-2, 2)) +
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
