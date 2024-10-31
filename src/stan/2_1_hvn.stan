data {
real si;
real b0; real bx; real by;
}
parameters {
  real x;
  real y;
 }
model {
real g = b0+bx*x+by*y;
real dgx = bx;real dgy = by;
real ndg =1;
target += normal_lpdf(0.00 | g/ndg, si); 
}
