function spatial_plot(G,theta,phi,bump_height,ref_sphere,plottype,edgeColor)
	switch plottype
		case 0
			F=abs(G).^2;
		case 1
			F=abs(G);
		case 2
			F=real(G);
		case 3
			F=imag(G);
		otherwise
			F=real(G);
	end
	maxF=max(max(abs(F)));
	F=F/maxF; % normalize entries to interval [-1.0,1.0]

	% make radius<=1 an affine mapping of the spherical harmonic value
	radius=abs(ref_sphere + bump_height*F)/(ref_sphere+bump_height);

	% convert to 3D Cartesian mesh
	rsint=radius.*sin(theta);
	x=rsint.*cos(phi);
	y=rsint.*sin(phi);
	z=radius.*cos(theta);

	% set up for cropping/saving to square image
	subplot('position',[-0.15 -0.15 1.3 1.3]);

	s=surf(x,y,z,F); % last argument determines colormap
	s.LineWidth=0.1; % for finer lines in printing/png
	s.EdgeColor=edgeColor; % no lines

	fig=gcf;
	fig.Position(3)=450;
	fig.Position(4)=450;

	fa=fig.CurrentAxes;

	light; lightangle(260,-45) % add 2 lights
	lighting gouraud % preferred lighting for a curved surface
	view(40,30) % set viewpoint

	maxa=1.0; % max of what radius above can be
	axis([-maxa maxa -maxa maxa -maxa maxa]);
	axis off

	c=colorbar('position',[0.01 0.75 0.03 0.22]);
	c.Ticks=[-1:0.5:1];
 	c.TickLabelInterpreter='latex';
	c.TickLabels={'$-1$','','$0$','','$1$'};
	c.FontSize=12;
 	c.AxisLocation='in';

	switch plottype
		case {0,1}
			c.Limits=[0,1];
		otherwise
			c.Limits=[-1,1];
	end

 	camzoom(1.15) % zoom into scene
 	axes(fa); rotate3d on % setup for for gui interaction

 	% output to png file to existing figures directory
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 6]) %150dpi
	saveas(gcf,'test','png')
end
