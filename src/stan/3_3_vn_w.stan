data {
  real<lower=0> si;
  real b0; real bx; real by; real bz;
  real bx2; real by2; real bz2;
  real bxy; real bxz; real byz;
  real bx3; real by3; real bz3;
  real bx2y; real bxy2; real bx2z; real bxz2; real by2z; real byz2; real bxyz;
  real w;
}

parameters {
  real<lower=-w,upper=w> x;
  real<lower=-w,upper=w> y;
  real<lower=-w,upper=w> z;
}

model {
  real g = b0 + bx*x + by*y + bx*z +  bx2*x^2 +  by2*y^2 + bz2*z^2+  bxy*x*y + bxz*x*z + byz*y*z + bx3*x^3 + by3*y^3 + bz3*z^3 + bx2y*x^2*y + bxy2*x*y^2 + bx2z*x^2*z + bxz2*x*z^2 +by2z*y^2*z + byz2*y*z^2 + bxyz*x*y*z;
  real ndg = 1;
  target += normal_lpdf(0.00 | g/ndg, si);
}
