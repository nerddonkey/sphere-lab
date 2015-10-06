function [s]=spatialPlot(G,theta,phi,bump_height,ref_sphere,plottype,doScale)
%spatialPlot

%% default parameters
if nargin<4
	bump_height=0.15;
end
if nargin<5
	ref_sphere=1.0;
end
if nargin<6
	plottype=1;
end
if nargin<7
	doScale=false;
end

F=spatialMap(G,plottype,doScale);

%% make radius<=1 an affine mapping of the spherical harmonic value
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

fig=gcf;
fig.Position(3)=600;
fig.Position(4)=600;

fa=fig.CurrentAxes;

delete(findall(gcf,'Type','light')) % remove existing lights
light; lightangle(260,-45) % add 2 lights
lighting gouraud % preferred lighting for a curved surface
view(40,30) % set viewpoint

maxa=1.0; % max of what radius above can be
axis([-maxa maxa -maxa maxa -maxa maxa]);
axis off

delete(findall(gcf,'Type','colorbar')) % remove existing colorbars
c=colorbar('position',[0.03 0.75 0.03 0.22]);
c.Ticks=-1:0.5:1;
c.TickLabelInterpreter='latex';
c.TickLabels={'$-1$','','$0$','','$1$'};
c.FontSize=12;
c.AxisLocation='in';

switch plottype
	case {0,1}
		c.Limits=[-0.02,1.02];
	otherwise
		c.Limits=[-1.02,1.02];
end

set(fa,'CameraViewAngle',9) % zoom into scene
axes(fa);
rotate3d on % setup for for gui interaction
