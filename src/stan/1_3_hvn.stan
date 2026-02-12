data {
  real si;
  real b1;   real bx;   real bx2;   real bx3;
}
parameters {
  real x;
}
transformed parameters {
  real g = b1+bx*x+bx2*x^2+bx3*x^3;
  real dgx = bx+2*bx2*x+3*bx3*x^2;
  real ndg = 1;
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
