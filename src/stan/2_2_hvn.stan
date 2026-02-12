data {
  real si;
  real b1;   real bx;   real by;   real bx2;   real by2;   real bxy;
}
parameters {
  real x;
  real y;
}
transformed parameters {
  real g = b1+bx*x+by*y+bx2*x^2+by2*y^2+bxy*x*y;
  real dgx = bx+2*bx2*x+bxy*y;
  real dgy = by+2*by2*y+bxy*x;
  real ndg = 1;
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
