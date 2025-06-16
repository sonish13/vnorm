data {
  real si;
  real b0;   real bx;   real by;
}
parameters {
  real x;
  real y;
}
transformed parameters {
  real g = b0+bx*x+by*y;
  real dgx = bx;
  real dgy = by;
  real ndg = 1;
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
