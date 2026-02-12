data {
  real si;
  real b1;   real bx;   real by;   real bz;
  real w;
}
parameters {
  real<lower=-w, upper=w> x;
  real<lower=-w, upper=w> y;
  real<lower=-w, upper=w> z;
}
transformed parameters {
  real g = b1+bx*x+by*y+bz*z;
  real dgx = bx;
  real dgy = by;
  real dgz = bz;
  real ndg = sqrt(dgx^2 + dgy^2 + dgz^2);
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
