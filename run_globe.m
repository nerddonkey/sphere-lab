function run_Globe
%run_Globe recontruct spherical data from data files found at
% https://geodesy.curtin.edu.au/research/models/Earth2012/
% http://www.ipgp.fr/~wieczor/SH/SH.html
% This demonstrates the inverse sphericsal harmonic transform, that is,
% going from the spectral domain to the spatial domain

close all;
L_max=600; % maximum included spherical harmonic degree
ntt=max(21,L_max+1); % number of points in theta
npp=max(41,2*L_max+1); % number of points in phi
% ntt=201;
% npp=221;

for bodyIndex=1:6
	switch bodyIndex
		case 1
			globe='Earth2012.topo_bathy_bed.SHCto2160.shape';
			name='Earth2012-bed';
			cmap='parula';
		case 2
			globe='Earth2012.topo_bathy.SHCto2160.shape';
			name='Earth2012';
			cmap='parula';
		case 3
			globe='Earth2012.topo_bathy.SHCto2160.shape';
			name='Earth2012_air';
			cmap='parula';
		case 4
			globe='MarsTopo2600.shape';
			name='Mars';
			cmap='hot';
		case 5
			globe='VenusTopo719.shape';
			name='Venus';
			cmap='hot';
		case 6
			globe='MoonTopo2600p.shape';
			name='Moon';
			cmap='parula';
	end

	[f,theta,phi]=spatial_globe(L_max,ntt,npp,globe);
	s=spatial_plot(f,theta,phi,0.05,1.0,2);
	colormap(cmap)
	s.EdgeColor='none'; % no lines

	llabel=sprintf('$L_{\\mathrm{max}}=%d$',L_max);
	delete(findall(gcf,'Tag','myLabel'));
	a=annotation('textbox',[0.7,0.89,0.28,0.1],'String',llabel);
	a.Interpreter='latex';
	a.FontSize=18;
	a.LineStyle='none';
	a.Tag='myLabel';

	fig=gcf;
	fig.Name=name;
	fig.Position(3)=600;
	fig.Position(4)=600;

	% output to png file to current directory
	outname=sprintf('figs/%s_%04d',name,L_max);
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
	saveas(gcf,outname,'png')
end

if 0 % render out frames for a movie on osx
	view(40,-15) % set viewpoint
	for i=0:359
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
		rotname=sprintf('frames/globerot_%04d', i);
		camorbit(1.0,0.0); drawnow; saveas(gcf,rotname,'png')
	end

	% setenv('DYLD_LIBRARY_PATH','/usr/local/bin/');
	system('/opt/local/bin/convert -delay 50 -loop 0 frames/globerot_*.png globerot.gif');
end
