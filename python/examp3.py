// example 3 : define then plot a spherical harmonic
// 3-1 : define the function Ylm
function [y]=Y(l, m, theta, phi)
  // theta may be a scalar or a row vector
  // phi may be a scalar or a column vector
  if m >= 0 then
     y = (-1)^m/(sqrt(2*%pi))*exp(%i*m*phi)*legendre(l, m, cos(theta), "norm")
  else
     y = 1/(sqrt(2*%pi))*exp(%i*m*phi)*legendre(l, -m, cos(theta), "norm")
  end
endfunction

// 3.2 : define another useful function
function [x, y, z]=sph2cart(theta, phi, r)
  // theta row vector      1 x nt
  // phi   column vector  np x 1
  // r     scalar or np x nt matrix (r(i,j) the length at phi(i) theta(j))
  x = r.*(cos(phi)*sin(theta));
  y = r.*(sin(phi)*sin(theta));
  z = r.*(ones(phi)*cos(theta));
endfunction

// 3-3 plot Y31(theta,phi)
l = 3; m = 1;
theta = linspace(0.1,%pi-0.1,60);
phi = linspace(0,2*%pi,120)';
f = Y(l,m,theta,phi);
[x1,y1,z1] = sph2cart(theta,phi,abs(f));       [xf1,yf1,zf1] = nf3d(x1,y1,z1);
[x2,y2,z2] = sph2cart(theta,phi,abs(real(f))); [xf2,yf2,zf2] = nf3d(x2,y2,z2);
[x3,y3,z3] = sph2cart(theta,phi,abs(imag(f))); [xf3,yf3,zf3] = nf3d(x3,y3,z3);

clf()
subplot(1,3,1)
plot3d(xf1,yf1,zf1,flag=[2 4 4]); xtitle("|Y31(theta,phi)|")
subplot(1,3,2)
plot3d(xf2,yf2,zf2,flag=[2 4 4]); xtitle("|Real(Y31(theta,phi))|")
subplot(1,3,3)
plot3d(xf3,yf3,zf3,flag=[2 4 4]); xtitle("|Imag(Y31(theta,phi))|")
