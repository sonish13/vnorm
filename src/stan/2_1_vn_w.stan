data {
  real si;
  real b1;   real bx;   real by;
  real w;
}
parameters {
  real<lower=-w, upper=w> x;
  real<lower=-w, upper=w> y;
}
transformed parameters {
  real g = b1+bx*x+by*y;
  real dgx = bx;
  real dgy = by;
  real ndg = sqrt(dgx^2 + dgy^2);
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
