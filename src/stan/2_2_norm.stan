data {
  real<lower=0> si;
  real b0; real bx; real by;
  real bx2; real bxy; real by2;
  real w;
}

parameters {
  real<lower=-w,upper=w> x;
  real<lower=-w,upper=w> y;
}

model {
  real g = b0 + bx*x + bx2*x^2 + bxy*x*y + by*y + by2*y^2;
  real dgx = bx + 2*bx2*x + bxy*y;
  real dgy = by + 2*by2*y + bxy*x;
  real ndg = sqrt(dgx^2 + dgy^2);
  target += normal_lpdf(0.00 | g/ndg, si);
}
