data {
  real si;
  real b1;   real bx;   real bx2;
  real w;
}
parameters {
  real<lower=-w, upper=w> x;
}
transformed parameters {
  real g = b1+bx*x+bx2*x^2;
  real dgx = bx+2*bx2*x;
  real ndg = 1;
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
