data {
  real si;
  real b1;   real bx;
}
parameters {
  real x;
}
transformed parameters {
  real g = b1+bx*x;
  real dgx = bx;
  real ndg = sqrt(dgx^2);
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
