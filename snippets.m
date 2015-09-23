
% sphereHarm:
% Y_l^m = (-1)^m Q_l P_l^m(sch) e^m\phi, for m\geq0
% Y_l^{-m} = (-1)^m conj(Y_l^m), for m>0
% trapSphereR: integration of region R defined by tt and pp
% \int_R Y_p^q(\unit{x}) \conj{Y_l^m(\unit{x})} ds(\unit{x})

% convert -delay 50 -loop 0 slep*.png slep.gif
% convert -delay 30 -quality 100 aus_* movie.mov
% check=bsxfun(@rdivide,D*w,w) % check eigenvector works as expected

setenv('PATH', [getenv('PATH') ';D:\Perl\bin']);

setenv('PATH', [getenv('PATH') ':/usr/local/bin:/usr/local/bin']);

Grid means it is expecting vectors tv, pv and forms a cartesian product grid theta,phi
by using [theta,phi]=ndgrid(tv,pv).

Spherical Harmonics on a Rectangular Grid
function [Ylm,theta,phi]=sphHarmGrid(l,m,tv,pv) computes Y_l^m
function [YB,theta,phi]=sphHarmGridBank(l,tv,pv) generates an cell array "Bank" Y_l^0,Y_l^1,...,Y_l^m
function [Ylm,theta,phi]=sphHarmGridBankm(l,m,tv,pv) same as sphHarmGrid but computed via sphHarmGridBank

Inverse Spherical Harmonic Transform


brew reinstall libav



%SPATIAL performs the inverse spherical harmonic transform. The spherical harmonic
% coefficients are provided in vector w which is organized in the (l,m): (0,0)
% (1,-1) (1,0) (1,1) (2,-2) ... ordering.  Inputs tv and pv specify the theta and
% phi values over which the spatial function f is computed.  The coefficient vector
% (either row or column) will be zero padded to be a squared dimension long.
%
% The first case tv and pv are both vectors (either row or column) which results in
% the computation of a separable mesh, or rectangular grid, using built-in ndgrid.
% Function f, being the complex-valued weighted sum of complex spherical harmonics,
% is computed over this theta-phi mesh.  As well as f, the created mesh theta and
% phi is returned. Computations are efficient and vectorized.
%
% The second case tv and pv are both matrices of identical size representing a
% general, not necessarily rectangular, grid.  Theses are taken as the mesh and
% copied to the output mesh theta and phi.  Function f is computed individually at
% each point on the mesh.  Computations are essentially all scalar because we can't
% assume separability.
%
% Y_l^m = (-1)^m Q_l P_l^m(sch) e^m\phi, for m\geq0
% Y_l^{-m} = (-1)^m conj(Y_l^m), for m>0
%
% SEPARABLE GRID CASE
% We exploit separability and the vector operations of matlab.
% The associated legendre function is done as a bank where rows
% are the order m and the columns are the vector of theta; this is
% for each degree l.  Similarly the complex exponentials are
% computed on a matrix mesh.  So the computation is pretty quick.

