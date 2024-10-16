data {
  real<lower=0> si;
  real b0; real bx; real by;
  real w;
}

parameters {
  real<lower=-w,upper=w> x;
  real<lower=-w,upper=w> y;
}

transformed parameters{
  real g = b0 + bx*x + by*y;
  real ndg = 1;
}

model {
  target += normal_lpdf(0.00 | g/ndg, si);
}
