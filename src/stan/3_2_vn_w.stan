data {
  real si;
  real b1;   real bx;   real by;   real bz;   real bx2;   real by2;   real bz2;   real bxy;   real bxz;   real byz;
  real w;
}
parameters {
  real<lower=-w, upper=w> x;
  real<lower=-w, upper=w> y;
  real<lower=-w, upper=w> z;
}
transformed parameters {
  real g = b1+bx*x+by*y+bz*z+bx2*x^2+by2*y^2+bz2*z^2+bxy*x*y+bxz*x*z+byz*y*z;
  real dgx = bx+2*bx2*x+bxy*y+bxz*z;
  real dgy = by+2*by2*y+bxy*x+byz*z;
  real dgz = bz+2*bz2*z+bxz*x+byz*y;
  real ndg = sqrt(dgx^2 + dgy^2 + dgz^2);
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
