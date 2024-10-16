data {
  real<lower=0> si;
  real b0; real bx; real by; real bz;
  real bx2; real by2; real bz2;
  real bxy; real bxz; real byz;
}

parameters {
  real x;
  real y;
  real z;
}

transformed parameters {
  real g = b0 + bx*x + by*y + bz * z + bx2*x^2 + by2*y^2 + bz2 * z^2 + bxy*x*y + bxz*x*z + byz*y*z;
  real ndg = 1;
}

model {
  target += normal_lpdf(0.00 | g/ndg, si);
}
