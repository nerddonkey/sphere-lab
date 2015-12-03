N=500;
tv = linspace(0, pi, N + 1);
pv = linspace(0, 2*pi, 2*N + 1);
tic;
w=0.5*ones(40000,1)*(1+1j);
[F,theta,phi]=ishtRectGrid(w,tv,pv,true);
toc;
