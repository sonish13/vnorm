data {
  real si;
  real b0;   real bx;   real by;   real bx2;   real by2;   real bxy;   real bx3;   real by3;   real bx2y;   real bxy2;
  real w;
}
parameters {
  real<lower=-w, upper=w> x;
  real<lower=-w, upper=w> y;
}
transformed parameters {
  real g = b0+bx*x+by*y+bx2*x^2+by2*y^2+bxy*x*y+bx3*x^3+by3*y^3+bx2y*x^2*y+bxy2*x*y^2;
  real dgx = bx+2*bx2*x+bxy*y+3*bx3*x^2+2*bx2y*x*y+bxy2*y^2;
  real dgy = by+2*by2*y+bxy*x+3*by3*y^2+bx2y*x^2+2*bxy2*x*y;
  real ndg = 1;
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
