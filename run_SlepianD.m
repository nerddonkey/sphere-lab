function run_SlepianD
%run_SlepianD compute and plot dominant Slepian eigenfunction on
% the sphere for a hardwired region and hardwired bandwidth
   L_max=15;
   deginc=2.0;

	% Sub-region 1 R2
   tt1=[30:deginc:60]*pi/180;
   pp1=[-70:deginc:-10]*pi/180;
 	D1=SlepianD(L_max,tt1,pp1);

	% Sub-region 2 R1
   tt2=[20:deginc:70]*pi/180;
   pp2=[0:deginc:40]*pi/180;
	D2=SlepianD(L_max,tt2,pp2);

	% D matrix for the union of (disjoint) subregions R=R1+R2
   D=D1+D2;
	fprintf('Size of D matrix: %d x %d\n', size(D))

	% dominant eigenvector in spectral domain
	[V,lambda]=eig(D);
	fprintf('Largest Eigenvalue: %8.6f\n', lambda(end,end))
	w=V(:,end);

	% dominant eigenfunction in spatial domain
	tt=[0:deginc:180]*pi/180;
	pp=[0:deginc:360]*pi/180;
	[slepian,theta,phi]=spatial(w,tt,pp,0);
	maxSlepian=max(abs(slepian(:)));
	slepian=slepian/maxSlepian; % normalize entries to interval [-1.0,1.0]

	% plot eigenfunction in spatial domain
	close
	plottype=1;
	bump_height=0.2;
	ref_sphere=1.0;
	s=spatial_plot(slepian,theta,phi,bump_height,ref_sphere,plottype);
	s.EdgeColor='none'; % no lines
	s.FaceAlpha=0.4;
	colormap(copper)

	hold on

	% highlight eigenfunction with region R
	[slepian1,theta,phi]=spatial(w,tt1,pp1,0);
	slepian1=slepian1/maxSlepian;
	spatial_plot(slepian1,theta,phi,bump_height,ref_sphere,plottype);

	[slepian2,theta,phi]=spatial(w,tt2,pp2,0);
	slepian2=slepian2/maxSlepian;
	spatial_plot(slepian2,theta,phi,bump_height,ref_sphere,plottype);

	llabel=sprintf('$L_{\\mathrm{max}}=%d$ $\\lambda=%8.6f$', L_max, lambda(end,end));
	delete(findall(gcf,'Tag','myLabel'));
	a=annotation('textbox',[0.71,0.89,0.28,0.1],'String',llabel);
	a.Interpreter='latex';
	a.FontSize=20;
	a.LineStyle='none';
	a.Tag='myLabel';

	hold off

 	% output to png file to current directory
	outname=sprintf('slep_%02d', L_max);
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
	saveas(gcf,outname,'png')
end

% convert -delay 50 -loop 0 slep*.png slep.gif
% check=bsxfun(@rdivide,D*w,w) % check eigenvector works as expected
