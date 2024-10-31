data {
real si;
real b0; real bx;
}
parameters {
  real x;
 }
model {
real g = b0+bx*x;
real dgx = bx;
real ndg = sqrt(dgx^2);
target += normal_lpdf(0.00 | g/ndg, si); 
}
