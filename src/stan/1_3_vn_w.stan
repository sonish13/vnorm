data {
  real<lower=0> si;
  real b0; real bx;
  real bx2; real bx3;
  real w;
}

parameters {
  real <lower=-w, upper = w>x;
}

transformed parameters{
  real g = b0 + bx *x + bx2 * x^2 + bx3 * x^3;
  real ndg = 1;
}

model {
  target += normal_lpdf(0.00 | g/ndg, si);
}
