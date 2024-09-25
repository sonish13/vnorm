data {
  real<lower=0> si;
  real b0; real bx; real by;
  real bx2; real bxy; real by2;
  real bx3; real bx2y; real bxy2; real by3;
  real w;
}

parameters {
  real<lower=-w,upper=w> x;
  real<lower=-w,upper=w> y;
}

model {
  real g = b0 + bx*x + bx2*x^2 + bxy*x*y + by*y + by2*y^2 + bx3 * x^3 + bx2y * x^2 * y + bxy2* x * y^2 + by3 * y^3;
  real dgx = bx + 2*bx2*x + bxy*y + 3 * bx3 * x^2 + 2 * bx2y * x * y + bxy2 * y^2;
  real dgy = by + 2*by2*y + bxy*x + 3 * by3 * y^2 + 2 * bxy2 * x * y + bx2y * x^2;
  real ndg = sqrt(dgx^2 + dgy^2);
  target += normal_lpdf(0.00 | g/ndg, si);
}
