function [F,theta,phi]=spatial_shn(n,tt,pp)
%SPATIAL_SHLM spatial recontruction of spherical harmonic of
% index n \in\{0,1,2,\dotsc,N_tot-1}

N_tot=ceil(sqrt(n+1))^2; % max number of coefficients
w=zeros(N_tot,1); % allocate for complex SH coefficients
w(n+1)=1;
[F,theta,phi]=spatial(w,tt,pp,0);
