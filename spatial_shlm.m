function [F,theta,phi]=spatial_shlm(l,m,tt,pp)
%SPATIAL_SHLM spatial recontruction of spherical harmonic of
% degree l and order m
	n=l*(l+1)+m; % (7.39)
	[F,theta,phi]=shn_spatial(n,tt,pp,0);
end
