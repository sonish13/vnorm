data {
  real si;
  real b0;   real bx;   real by;   real bz;
}
parameters {
  real x;
  real y;
  real z;
}
transformed parameters {
  real g = b0+bx*x+by*y+bz*z;
  real dgx = bx;
  real dgy = by;
  real dgz = bz;
  real ndg = 1;
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
