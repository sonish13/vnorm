data {
  real si;
  real b0;   real bx;   real by;   real bx2;   real by2;   real bxy;
  real w;
}
parameters {
  real<lower=-w, upper=w> x;
  real<lower=-w, upper=w> y;
}
transformed parameters {
  real g = b0+bx*x+by*y+bx2*x^2+by2*y^2+bxy*x*y;
  real dgx = bx+2*bx2*x+bxy*y;
  real dgy = by+2*by2*y+bxy*x;
  real ndg = sqrt(dgx^2 + dgy^2);
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
