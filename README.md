# sphere-lab

## Matlab code for spherical harmonic computations

Some routines for spherical harmonic transform related work

### Spherical Harmonic

- function [Ylm,theta,phi]=sphereHarm(l,m,tt,pp)
- function [I]=trapSphereR(f,tt,pp) &mdash; \int_R f(x) ds(x) using 2D trapezoid rule integration in two lines
- function [f,theta,phi]=spatial(w,tt,pp,useProgressBar) &mdash; performs the inverse spherical harmonic transform
- function [y,theta,phi]=SlepianDrc(r,c,tt,pp) &mdash; computes D_{r,c} in (8.27) where r={l,m} (row) and c={p,q} (column)
