
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
and “standard deviation” equal to .

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
#> Total execution time: 1.3 seconds.
head(samps)
#>            x         y
#> 1 -0.8139800  0.415087
#> 2 -0.0591321  0.979381
#> 3  0.9203630 -0.294203
#> 4  0.9970770 -0.281915
#> 5  0.1841940  1.022640
#> 6  0.1908310  0.923108
str(samps)
#> 'data.frame':    2000 obs. of  2 variables:
#>  $ x: num  -0.814 -0.0591 0.9204 0.9971 0.1842 ...
#>  $ y: num  0.415 0.979 -0.294 -0.282 1.023 ...
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
#> Chain 1 finished in 0.4 seconds.
#> Chain 2 finished in 0.3 seconds.
#> Chain 3 finished in 0.3 seconds.
#> Chain 4 finished in 0.3 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.3 seconds.
#> Total execution time: 1.5 seconds.
head(samps)
#>           x        y          z
#> 1  0.453262 0.799288 -0.0800479
#> 2  0.443706 0.803808 -0.0160269
#> 3 -0.239365 0.827039  0.3014040
#> 4 -0.847219 0.303297  0.6017080
#> 5 -0.776019 0.300545  0.5745230
#> 6 -0.764296 0.306352  0.5701080
```

We also want the user to be able to pre-compile models. This can be done
with `compile_stan_code`. This helps avoid re-compiling stan model for
similar polynomial with different coefficients.

``` r
poly <- mpoly::mp("x^4 + y^4 - 1")
compile_stan_code(poly = poly)
#> compiled_stan_info variable created in Global environment
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
#> Chain 1 finished in 0.0 seconds.
#> Chain 2 finished in 0.1 seconds.
#> Chain 3 finished in 0.0 seconds.
#> Chain 4 finished in 0.0 seconds.
#> 
#> All 4 chains finished successfully.
#> Mean chain execution time: 0.0 seconds.
#> Total execution time: 0.5 seconds.
head(samps)
#>           x         y
#> 1 -0.284020  1.135140
#> 2 -0.385049  0.841975
#> 3 -0.856798 -0.536129
#> 4 -0.893442 -0.255433
#> 5 -0.801538  0.354110
#> 6  0.180647 -0.730491
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
