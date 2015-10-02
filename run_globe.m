function run_Globe
%%run_Globe recontruct spherical data from shape planetary data files
% found at:
%    https://geodesy.curtin.edu.au/research/models/Earth2012/
%    http://www.ipgp.fr/~wieczor/SH/SH.html
%
% This demonstrates the inverse sphericsal harmonic transform, that is,
% going from the spectral domain to the spatial domain.

%%
L_max=300; % maximum included spherical harmonic degree
ntt=max(21,L_max+1); % number of points in theta
npp=max(41,2*L_max+1); % number of points in phi (here first and last phi are the same)
bump=0.05;
doScale=1;
radius=0.0; % use file value

%% Iterate over celestial bodies (files need to be downloaded and in the path)
for bodyIndex=2:2
	switch bodyIndex
		case 1
			globe='Earth2012.topo_bathy_bed.SHCto2160.shape';
			name='Earth2012-bed';
			cmap='parula';
			radius=6371000;
			bump=0.035;
			az=40;  el=-15;
		case 2
			globe='Earth2012.topo_bathy.SHCto2160.shape';
			name='Earth2012';
			cmap='parula';
			radius=6371000;
			bump=0.035;
			az=40;  el=-15;
		case 3
			globe='Earth2012.topo_air.SHCto2160.shape';
			name='Earth2012_air';
			cmap='parula';
			radius=6371000;
			bump=0.035;
			az=40;  el=-15;
		case 4
			globe='MarsTopo2600.shape';
			name='Mars';
			cmap='parula';
			bump=0.05;
			az=180;  el=0;
		case 5
			globe='VenusTopo719.shape';
			name='Venus';
			cmap='parula';%'hot';
			bump=0.035;
			az=180;  el=0;
		case 6
			globe='MoonTopo2600p.shape';
			name='Moon';
			cmap='parula';
			bump=0.02;
			az=180;  el=0;
	end

	%% Get the data from the shape file
	[F,theta,phi,L_max,rad]=ishtFromShapeFile(L_max,ntt,npp,globe,doScale);

	if radius==0
		radius=rad;
	end

	% figure and data output
	base=userpath;  base(end)='/'; % ~/Documents/MATLAB
	figs_folder=[base 'figs/'];
	frames_basename=sprintf('%s_%04d',[base 'frames/' name],L_max);
	figs_basename=sprintf('%s_%04d',[figs_folder name],L_max);
	save([figs_basename '.mat'],'F','theta','phi');

	% thinking of making the rest a function so the F,theta,phi data in a *.mat
	% file can be drawn
	% renderGlobe(F,theta,phi,L_max,name,0.05,1.0,2,cmap,az,el,figs_basename,frames_basename)

	%% Render the globe to figure
	close all;
	s=spatialPlot(F,theta,phi,bump,1.0,2);%0.05
	colormap(cmap)
	s.EdgeColor='none'; % no lines

	%% Annotate the figure
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

	%% Output png file to current directory
	set(gcf,'InvertHardCopy','off');
	set(gcf,'color',[0.7 0.7 0.7]); % Set the figure frame color to white
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
	saveas(gcf,figs_basename,'png')

	%% Render movie on osx; needs avconv - get via 'brew install libav'
	if ~system('which avconv >/dev/null')
		view(az,el) % set viewpoint
		for i=0:358
			set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
			rotname=sprintf('%s_%04d',frames_basename,i);
			camorbit(1.0,0.0); drawnow; saveas(gcf,rotname,'png')
		end
		renderMovWithAvconv=sprintf(...
			'avconv -framerate 30 -y -v quiet -f image2 -i %s_%%04d.png %s.mov',...
			frames_basename,frames_basename);
		if ~system(renderMovWithAvconv) % make compressed mov
			cleanup=sprintf('rm %s*.png',frames_basename);
			system(cleanup); % delete frames
		end
	end
end
