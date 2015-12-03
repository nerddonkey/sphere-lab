N=900;
l=25;
m=2;
tv = linspace(0, pi, N + 1);
pv = linspace(0, 2*pi, 2*N + 1);
[theta,phi]=ndgrid(tv,pv);
tic;
Ylm=sphHarm(l,m,theta,phi);
toc;

[l,m,size(Ylm)]