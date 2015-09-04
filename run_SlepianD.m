function run_SlepianD
%run_SlepianD compute and plot dominant Slepian eigenfunction on
% the sphere for a hardwired region and hardwired bandwidth
for eigindex=0:0 % eigenvalue index 0,1,2 in descending energy order
for L_max=6:6 % range of L_max
	intginc=0.1; % fine grid for integration
	medinc=0.5; % coarse grid for plotting
	plotinc=2.0; % coarse grid for plotting

   fprintf('\n@@ L_max: %d\n',L_max)
   fprintf('@@ Fine Grid: %.2f\n',intginc)
   fprintf('@@ Coarse Grid: %.2f\n',plotinc)

	% (disjoint) subregions (3 off)
% 	tr=[30 60; 20 70; 80 110];
% 	pr=[-70 -10; 0 40; -70 -10];
	tr=[30 60; 80 110];
	pr=[-70 -10; -70 -10];

	N_max=(L_max+1)^2;
	D=zeros(N_max,N_max);

	for r=1:size(tr,1) % loop over subregions
		tt=(tr(r,1):intginc:tr(r,2))*pi/180;
		pp=(pr(r,1):intginc:pr(r,2))*pi/180;
		D=D+SlepianDH(L_max,tt,pp);
      fprintf('@@ Completed D for Region: %d\n',r)
	end

	fprintf('@@ Size of D matrix: %dx%d\n', size(D))

	% dominant eigenvector w with eigenvalue lambda in spectral domain
	[V,lamD]=eig(D); % eigenvalue
	lambda=lamD(end-eigindex,end-eigindex);
   % stem(flip(diag(lamD))) % plot the eigenvalues
	fprintf('@@ Largest Eigenvalue: %8.6f\n',lambda)
	w=V(:,end-eigindex); % eigenvector
	w=w/norm(w); % normalize eigenvector (superfluous)

	% dominant eigenfunction in spatial domain
	tt=(0:intginc:180)*pi/180;
	pp=(0:intginc:360)*pi/180;
	[slepian,theta,phi]=spatial(w,tt,pp,0);
	maxSlepian=1.2;%max(abs(slepian(:))); % normalization to [-1.0,1.0] for plotting
	fprintf('@@ Max Slepian: %8.6f\n',maxSlepian)

	% energies in spatial and spectral domains
	spectralEnergy=norm(w)^2;
	spatialEnergy=trapSphereR(abs(slepian).^2,tt,pp);
	fprintf('@@ Spectral Energy: %8.6f\n',spectralEnergy)
	fprintf('@@ Spatial Energy: %8.6f\n',spatialEnergy)

	% plot eigenfunction in spatial domain
	close
	bump_height=0.15;
	ref_sphere=1.0;
	plottype=1;

 	tt=(0:medinc:180)*pi/180;
 	pp=(0:medinc:360)*pi/180;
 	[slepian,theta,phi]=spatial(w,tt,pp,0);
	s=spatial_plot(slepian/maxSlepian,theta,phi,bump_height,ref_sphere,plottype);
	s.EdgeColor='none'; % no lines
	s.FaceAlpha=0.8;
	colormap(copper)

	hold on

	% Compute energy in sub-regions
	lambda_est=0;
	for r=1:size(tr,1) % loop over subregions (fine grid)
		tt=(tr(r,1):intginc:tr(r,2))*pi/180;
		pp=(pr(r,1):intginc:pr(r,2))*pi/180;
		reg=spatial(w,tt,pp,0);
		lambda_est=lambda_est+trapSphereR(abs(reg).^2,tt,pp);
   end
   lambda_est=lambda_est/spatialEnergy;

	% Plot sub-regions computation
	for r=1:size(tr,1) % loop over subregions (coarse grid)
		tt=(tr(r,1):plotinc:tr(r,2))*pi/180;
		pp=(pr(r,1):plotinc:pr(r,2))*pi/180;
		[reg,theta,phi]=spatial(w,tt,pp,0);
		spatial_plot(reg/maxSlepian,theta,phi,bump_height,ref_sphere,plottype);
	end

	llabel=sprintf('$L_{\\mathrm{max}}=%d$ $\\lambda=%8.6f$',L_max,lambda);
	delete(findall(gcf,'Tag','myLabel'));
	a=annotation('textbox',[0.71,0.89,0.28,0.1],'String',llabel);
	a.Interpreter='latex';
	a.FontSize=20;
	a.LineStyle='none';
	a.Tag='myLabel';

	fprintf('@@ Region Energy: %8.6f (eigenvalue %8.6f)\n',lambda_est,lambda)

	hold off

	% output to png file (check directory exists)
	outname=sprintf('figs/sleppy_%02d_%02d',L_max,eigindex);
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
	saveas(gcf,outname,'png')
end
end
end

% convert -delay 50 -loop 0 slep*.png slep.gif
% check=bsxfun(@rdivide,D*w,w) % check eigenvector works as expected
