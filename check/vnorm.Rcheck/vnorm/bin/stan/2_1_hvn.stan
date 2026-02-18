data {
  real si;
  real b1;   real bx;   real by;
}
parameters {
  real x;
  real y;
}
transformed parameters {
  real g = b1+bx*x+by*y;
  real dgx = bx;
  real dgy = by;
  real ndg = 1;
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
