data {
  real si;
  real b0;   real bx;   real by;   real bz;   real bx2;   real by2;   real bz2;   real bxy;   real bxz;   real byz;   real bx3;   real by3;   real bz3;   real bx2y;   real bxy2;   real bx2z;   real by2z;   real bxz2;   real byz2;   real bxyz;
}
parameters {
  real x;
  real y;
  real z;
}
transformed parameters {
  real g = b0+bx*x+by*y+bz*z+bx2*x^2+by2*y^2+bz2*z^2+bxy*x*y+bxz*x*z+byz*y*z+bx3*x^3+by3*y^3+bz3*z^3+bx2y*x^2*y+bxy2*x*y^2+bx2z*x^2*z+by2z*y^2*z+bxz2*x*z^2+byz2*y*z^2+bxyz*x*y*z;
  real dgx = bx+2*bx2*x+bxy*y+bxz*z+3*bx3*x^2+2*bx2y*x*y+bxy2*y^2+2*bx2z*x*z+bxz2*z^2+bxyz*y*z;
  real dgy = by+2*by2*y+bxy*x+byz*z+3*by3*y^2+bx2y*x^2+2*bxy2*x*y+2*by2z*y*z+byz2*z^2+bxyz*x*z;
  real dgz = bz+2*bz2*z+bxz*x+byz*y+3*bz3*z^2+bx2z*x^2+by2z*y^2+2*bxz2*x*z+2*byz2*y*z+bxyz*x*y;
  real ndg = sqrt(dgx^2 + dgy^2 + dgz^2);
}
model {
  target += normal_lpdf(0.00 | g/ndg, si); 
}
