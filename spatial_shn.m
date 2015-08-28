function [f,theta,phi]=spatial_shn(n,tt,pp)
%SPATIAL_SHLM spatial recontruction of spherical harmonic of
% index n \in\{0,1,2,\dotsc,N_max-1}
	N_max=ceil(sqrt(n))^2; % max number o
	w=zeros(N_max,1); % allocate for complex SH coefficients
	w(n)=1;
	[f,theta,phi]=spatial(w,tt,pp);
end
