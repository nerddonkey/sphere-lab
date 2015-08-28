function [y,f,theta,phi]=herm(n1,n2,tt,pp)
	[f,theta,phi]=spatial_shn(n1,tt,pp);
	[g,~,~]=spatial_shn(n2,tt,pp);
	f=f.*conj(g);
	y=sum(f(:));
% 	plot the mesh
% 	spatial_plot(f,theta,phi,0.5,1.0,2,'none');
% 	colormap('parula')
end
