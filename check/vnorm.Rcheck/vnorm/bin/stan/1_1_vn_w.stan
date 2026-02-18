data {
  real si;
  real b1;   real bx;
  real w;
}
parameters {
  real<lower=-w, upper=w> x;
}
transformed parameters {
  real g = b1+bx*x;
  real dgx = bx;
  real ndg = sqrt(dgx^2);
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
