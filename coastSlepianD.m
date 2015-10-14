function coastSlepianD(r,orthonorm,testEnergy)
%run_SlepianD compute and plot Slepian eigenfunctions on the sphere (8.27) and (8.29)

if nargin<1
	r=1; % default is to do Australia including Tasmania
end
if nargin<2
	orthonorm=false;
end
if nargin<3
	testEnergy=false;
end

%% Close all figures
close all

%% Coastline data from function regionIndexRange(minCoast)
coast(1).name='australia'; coast(1).range=[8296:8604 NaN 8623:8645]; coast(1).voff=[10 30];
coast(2).name='mainland';  coast(2).range=8296:8604; coast(2).voff=[10 30];
coast(3).name='tasmania';  coast(3).range=8623:8645; coast(3).voff=[10 30];
coast(4).name='eurasia';   coast(4).range=4796:6667; coast(4).voff=[0 0];
coast(5).name='antartica'; coast(5).range=0001:0612; coast(5).voff=[0 0];
coast(6).name='america';   coast(6).range=0741:2386; coast(6).voff=[0 0];
coast(7).name='africa';    coast(7).range=4281:4794; coast(7).voff=[0 0];
coast(8).name='greenland'; coast(8).range=3977:4210; coast(8).voff=[-20 20];
coast(9).name='newzealand'; coast(9).range=[8647:8738 NaN 8740:8786]; coast(9).voff=[10 20];
coast(10).name='anzac'; coast(10).range=[8296:8604 NaN 8623:8645 NaN 8647:8738 NaN 8740:8786]; coast(10).voff=[10 10];
if r>10; error('Invalid region'); end;

%% Determine paths for output files
base=userpath; base(end)='/'; % ~/Documents/MATLAB/
frames_folder=[base 'cslepian-frames/']; % ~/Documents/MATLAB/cslepian-frames/
dataFolder=[base 'cslepian-data/']; % ~/Documents/MATLAB/cslepian-data/
basename=coast(r).name;
diary(sprintf('%s_diary.txt',[dataFolder basename]))
diary on
disp(['@@ ',basename,' Started: ',datestr(now)])

thetaDeg=1.0;
phiDeg=1.0;
[thetaVecR,phiVecR,maskR,thetaR,phiR]=coastRegion(coast(r).range,thetaDeg,phiDeg);

tv=thetaVecR*pi/180;
pv=phiVecR*pi/180;
tR=thetaR*pi/180;
pR=phiR*pi/180;

%% Populate the Ntot x Ntot Hermitian D matrix (8.27)
for Lmax=20%:10:50
	fprintf('\n@@ Lmax: %d\n',Lmax)
	Ntot=(Lmax+1)^2;
	D=SlepianDH(Lmax,tv,pv,maskR,true);
	fprintf('@@ Size of D matrix: %dx%d\n', size(D))

	%% Save the D matrix
	dataName=sprintf('%s_%04d',[dataFolder basename],Lmax);
	save([dataName '.mat'],'Lmax','D');

	%% Get spectral eigenstructure and sort
	[V,lamD]=eig(D); % eigen-structure
	[lambda,order]=sort(diag(lamD),'descend'); % sort eigenvalues in descending order
	lambda=min(max(lambda,0.0),1.0); % pin to [0,1] (should be superfluous)
	V=V(:,order);

	stem(lambda) % plot the eigenvalues
	set(gca,'yscal','log'); shg; pause(2)

	maxSlepian=-1.0;

	for n=1:3%Ntot-1 % eigenvalue index in descending energy order
		if n>Ntot; break; end
		fprintf('@@ Eigenvalue %04d: %8.6f\n',n-1,lambda(n))

		%% n'th eigenvector w with eigenvalue lambda in spectral domain
		w=V(:,n); % eigenvector - SHT of Slepian function
		w=w/norm(w); % normalize eigenvector (actually this is superfluous)

		%% Plot eigenfunction in spatial domain
		close % prepare fresh figure

		%% Do the inverse SHT
		plotinc=0.2;
		tvPlot=(0:plotinc:180)*pi/180;
		pvPlot=(0:plotinc:360)*pi/180;
		pvPlot=pvPlot+pv(1); % (1) move phi seam to left hand edge so no discontinuity

		[slepianPlot,thetaPlot,phiPlot]=ishtRectGrid(w,tvPlot,pvPlot,true);

		if orthonorm
			slepianPlot=slepianPlot/sqrt(lambda(n));

			%% Plot only orthonormal functions in the region
	%		maskRPlot=inpolygon(phiPlot,thetaPlot,phiR*pi/180,thetaR*pi/180);
	%		slepianPlot=slepianPlot.*maskRPlot/sqrt(lambda(n));
		end

		if maxSlepian<0
			maxSlepian=1.0*max(abs(slepianPlot(:))); % for plot normalization
			fprintf('@@ Max Slepian: %8.6f\n',maxSlepian)
		end

		%% Plot the result
		bump=0.3; ref=1.0; ptype=1;
		s=spatialPlot(slepianPlot/maxSlepian,thetaPlot,phiPlot,bump,ref,ptype);

		%% Tweak the appearance
		s.EdgeColor='none'; % no lines
		s.FaceAlpha=0.8;
		colormap('cool') % or copper
		set(gcf,'MenuBar','none');
		set(gcf,'ToolBar','none');
		hold on

		%% Draw coastline on Slepian surface
		rR=interpn(thetaPlot,phiPlot,slepianPlot,tR,pR); % needs the rotation in (1)
		F=spatialMap(rR/maxSlepian,ptype);
		rad=abs(ref + bump*F)/(ref+bump)*1.02;
		coastR=[rad.*sin(tR).*cos(pR); rad.*sin(tR).*sin(pR); rad.*cos(tR)];
		plot3(coastR(1,:),coastR(2,:),coastR(3,:),'w','LineWidth',2);

		[az,el]=view(nanmean(coastR,2)); % centre of coastline
		view(az+coast(r).voff(2),el+coast(r).voff(1)); % offset for aesthetics
		set(gca,'CameraViewAngle',9) % zoom into scene 9

		%% Annotate plot
		llabel=sprintf('$L_{\\mathrm{max}}=%d$\n$\\lambda_{%d}=%8.6f$', Lmax,n-1,lambda(n));
		delete(findall(gcf,'Tag','myLabel'));
		a=annotation('textbox',[0.695,0.895,0.28,0.1],'String',llabel);
		a.Interpreter='latex';
		a.FontSize=18;
		a.LineStyle='none';
		a.Tag='myLabel';
		hold off

		%% Output to png file (directory needs to exist)
		set(gcf,'InvertHardCopy','off');
		set(gcf,'color','w'); % Set the figure frame color to white
		set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
		outname=sprintf('%s_%04d_%04d',[frames_folder basename],Lmax,n);
		saveas(gcf,outname,'png')
		if orthonorm && lambda(n)<0.001 break; end;
	end

	%% Render eigenfunction movie on osx; needs avconv - get via brew install libav
	if ~system('which avconv >/dev/null')
		renderMovWithAvconv=sprintf(...
			'avconv -framerate 2 -y -v quiet -f image2 -i %s_%04d_%%04d.png %s-%04d.mov',...
			[frames_folder basename],Lmax,[frames_folder basename],Lmax);
		system(renderMovWithAvconv); % make compressed mov
		cleanup=sprintf('rm %s_%04d*.png',[frames_folder basename],Lmax);
		system(cleanup); % delete frames
	end
end

disp(['@@ Ended: ' datestr(now)])
diary off
