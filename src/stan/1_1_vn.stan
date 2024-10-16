data {
  real<lower=0> si;
  real b0; real bx;
}

parameters {
  real x;
}

transformed parameters{
  real g = b0 + bx*x;
  real ndg = 1;
}

model {
  target += normal_lpdf(0.00 | g/ndg, si);
}
