function run_SlepianD
%run_SlepianD compute and plot Slepian eigenfunctions on the sphere (8.27) and (8.29)

intginc=0.1; % fine grid for integration (degrees)
medinc=0.5; % medium grid for smooth plotting (degrees)
plotinc=2.0; % coarse grid for grid plotting (degrees)

fprintf('\n@@ Fine Grid: %.2f (degrees)\n',intginc)
fprintf('@@ Medium Grid: %.2f (degrees)\n',medinc)
fprintf('@@ Coarse Grid: %.2f (degrees)\n',plotinc)

ttf=(0:intginc:180)*pi/180;
ppf=(0:intginc:360)*pi/180;

% subregion 1
tr(1,:)=[30 60];
pr(1,:)=[-70 -10];
% subregion 2
tr(2,:)=[80 110];
pr(2,:)=[-70 -10];
% subregion 3
tr(3,:)=[20 70];
pr(3,:)=[0 40];

for eigindex=0:0 % eigenvalue index 0,1,2 in descending energy order
	for L_max=4:4 % range of L_max
		fprintf('\n@@ L_max: %d\n',L_max)
		N_tot=(L_max+1)^2;
		D=zeros(N_tot,N_tot); % allocate and initialize to zero

		for r=1:size(tr,1) % loop over subregions
			ttfr=(tr(r,1):intginc:tr(r,2))*pi/180;
			ppfr=(pr(r,1):intginc:pr(r,2))*pi/180;
			D=D+SlepianDH(L_max,ttfr,ppfr);
			fprintf('@@ Completed D for Region: %d\n',r)
		end
		fprintf('@@ Size of D matrix: %dx%d\n', size(D))

		% dominant eigenvector w with eigenvalue lambda in spectral domain
		[V,lamD]=eig(D); % eigenvalue
		lambda=lamD(end-eigindex,end-eigindex);
		fprintf('@@ Eigenvalue %d: %8.6f\n',eigindex,lambda)
		w=V(:,end-eigindex); % eigenvector
		w=w/norm(w); % normalize eigenvector (superfluous)
		% stem(flip(diag(lamD))) % plot the eigenvalues

		% dominant eigenfunction in spatial domain
		slepian=spatial(w,ttf,ppf);
		maxSlepian=1.2;%max(abs(slepian(:))); % normalization to [-1.0,1.0] for plotting
		fprintf('@@ Max Slepian: %8.6f\n',maxSlepian)

		% energies in spatial and spectral domains
		spectralEnergy=norm(w)^2;
		spatialEnergy=trapSphereR(abs(slepian).^2,ttf,ppf);
		fprintf('@@ Spectral Energy: %8.6f\n',spectralEnergy)
		fprintf('@@ Spatial Energy: %8.6f\n',spatialEnergy)

		% plot eigenfunction in spatial domain
		close
		bump_height=0.15; ref_sphere=1.0; plottype=1;
		[slepian,theta,phi]=spatial(w,(0:medinc:180)*pi/180,(0:medinc:360)*pi/180,1);
		s=spatial_plot(slepian/maxSlepian,theta,phi,bump_height,ref_sphere,plottype);
		s.EdgeColor='none'; % no lines
		s.FaceAlpha=0.8;
		colormap(copper)
		hold on

		% Compute energy in sub-regions
		lambda_est=0;
		for r=1:size(tr,1) % loop over subregions (fine grid)
			ttfr=(tr(r,1):intginc:tr(r,2))*pi/180;
			ppfr=(pr(r,1):intginc:pr(r,2))*pi/180;
			reg=spatial(w,ttfr,ppfr,0);
			lambda_est=lambda_est+trapSphereR(abs(reg).^2,ttfr,ppfr);
		end
		lambda_est=lambda_est/spatialEnergy;
		fprintf('@@ Region Energy: %8.6f (eigenvalue %8.6f)\n',lambda_est,lambda)

		% Plot sub-regions
		for r=1:size(tr,1) % loop over subregions (coarse grid)
			ttcr=(tr(r,1):plotinc:tr(r,2))*pi/180;
			ppcr=(pr(r,1):plotinc:pr(r,2))*pi/180;
			[reg,theta,phi]=spatial(w,ttcr,ppcr,0);
			spatial_plot(reg/maxSlepian,theta,phi,bump_height,ref_sphere,plottype);
		end

		% Annotate plot
		llabel=sprintf('$L_{\\mathrm{max}}=%d$ $\\lambda=%8.6f$',L_max,lambda);
		delete(findall(gcf,'Tag','myLabel'));
		a=annotation('textbox',[0.71,0.89,0.28,0.1],'String',llabel);
		a.Interpreter='latex';
		a.FontSize=20;
		a.LineStyle='none';
		a.Tag='myLabel';
		hold off

		% output to png file (check directory exists)
		outname=sprintf('figs/sleppy_%02d_%02d',L_max,eigindex);
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
		saveas(gcf,outname,'png')
	end
end

% convert -delay 50 -loop 0 slep*.png slep.gif
% check=bsxfun(@rdivide,D*w,w) % check eigenvector works as expected
