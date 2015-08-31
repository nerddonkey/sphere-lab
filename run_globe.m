function []=run_Globe
%run_Globe
% https://geodesy.curtin.edu.au/research/models/Earth2012/
% http://www.ipgp.fr/~wieczor/SH/SH.html
	L_max=80  ; % maximum included spherical harmonic degree
	ntt=max(21,L_max+1); % number of points in theta
	npp=max(41,2*L_max+1); % number of points in phi
	%ntt=46;
	%npp=91;

% 	make sure data files directory is in path
	earth1='Earth2012.topo_bathy_bed.SHCto2160.shape';
	earth2='Earth2012.topo_bathy.SHCto2160.shape';
	earth3='Earth2012.topo_air.SHCto2160.shape';
	[f,theta,phi]=spatial_globe(L_max,ntt,npp,earth2);
	colormap('parula')
	s=spatial_plot(f,theta,phi,0.05,1.0,2);
	s.EdgeColor='none'; % no lines

% 	mars='MarsTopo2600.shape';
% 	[f,theta,phi]=spatial_globe(L_max,ntt,npp,mars);
% 	colormap('pink')
% 	spatial_plot(f,theta,phi,0.1,1.0,2,'none');

	venus='VenusTopo719.shape';

% 	moon='MoonTopo2600p.shape';
% 	[f,theta,phi]=spatial_globe(L_max,ntt,npp,moon);
% 	colormap('bone')
% 	spatial_plot(f,theta,phi,0.02,1.0,2,'none');

	llabel=sprintf('$L_{\\mathrm{max}}=%d$',L_max);
	delete(findall(gcf,'Tag','myLabel'));
	a=annotation('textbox',[0.71,0.89,0.28,0.1],'String',llabel);
	a.Interpreter='latex';
	a.FontSize=20;
	a.LineStyle='none';
	a.Tag='myLabel';

 	% output to png file to current directory
	outname=sprintf('figs/globe_%04d', L_max);
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
	saveas(gcf,outname,'png')

	%for i=[1:400]; camorbit(0.9,-0.1); drawnow; end
end
