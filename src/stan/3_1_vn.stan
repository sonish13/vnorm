data {
real si;
real b0; real bx; real by; real bz;
}
parameters {
  real x;
  real y;
  real z;
 }
model {
real g = b0+bx*x+by*y+bz*z;
real dgx = bx;real dgy = by;real dgz = bz;
real ndg = sqrt(dgx^2 + dgy^2 + dgz^2);
target += normal_lpdf(0.00 | g/ndg, si); 
}
