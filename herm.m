function [y,f,theta,phi]=herm(n1,n2,tt,pp)
	[f,theta,phi]=spatial_shn(n1,tt,pp);
	[g,~,~]=spatial_shn(n2,tt,pp);
	% integrand
	f=f.*conj(g);
	% apply the measure distortion
	f=bsxfun(@times,f,sin(tt(:)));
	% sum to approximate integration
	dttdpp=(tt(2)-tt(1))*(pp(2)-pp(1)); % uniform grid
	y=sum(f(:))*dttdpp;
% 	plot the mesh
% 	spatial_plot(f,theta,phi,0.5,1.0,2,'none');
% 	colormap('par  ula')
end
