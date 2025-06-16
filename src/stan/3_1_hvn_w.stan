data {
  real si;
  real b0;   real bx;   real by;   real bz;
  real w;
}
parameters {
  real<lower=-w, upper=w> x;
  real<lower=-w, upper=w> y;
  real<lower=-w, upper=w> z;
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
