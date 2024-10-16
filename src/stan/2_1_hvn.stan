data {
  real<lower=0> si;
  real b0; real bx; real by;
}

parameters {
  real x;
  real y;
}

transformed parameters{
  real g = b0 + bx*x + by*y;
  real dgx = bx;
  real dgy = by;
  real ndg = sqrt(dgx^2 + dgy^2);
}

model {
  target += normal_lpdf(0.00 | g/ndg, si);
}
