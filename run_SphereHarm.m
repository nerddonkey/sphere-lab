function run_SphereHarm
%run_SphereHarm tests sphereHarm

l=5;
m=1;
deginc=0.25;
tt=[0:deginc:180]*pi/180;
pp=[0:deginc:360]*pi/180;

[Y,theta,phi]=sphereHarm(l,m,tt,pp);

% normalize entries to interval [-1.0,1.0]
maxY=max(abs(Y(:)));
Y=Y/maxY;

bump_height=0.8; ref_sphere=1.0; plottype=2; % real
s=spatial_plot(Y,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='none'; % no lines
