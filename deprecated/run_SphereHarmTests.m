function run_SphereHarmTests
%run_SphereHarm tests various configurations of teh computation of
% spherical harmonics on separable and non-separable meshes

close all

l=5;
m=3;

% pick delta spherical harmonic coefficients for ISHT
n=l*(l+1)+m; % (7.39)
N_tot=ceil(sqrt(n+1))^2; % max number of coefficients
w=zeros(N_tot,1); % allocate for complex SH coefficients
w(n+1)=1;

deginc=2.0;
tv=(0:deginc:180)*pi/180;
pv=(0:deginc:360)*pi/180;

pos=50;pinc=200;


%%%% TEST separable sphereHarm

[Ylm,theta,phi]=sphereHarm(l,m,tv,pv);

doPlot(Ylm,theta,phi,pos,'separable sphereHarm');


%%%% TEST separable sphereHarmBankm/sphereHarmBank

[Ylm,theta,phi]=sphereHarmBankm(l,m,tv,pv);

pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'separable sphereHarmBankm');


%%%% TEST separable spatial (ISHT)

[Ylm,~,~]=spatial(w,tv,pv,1); % non-separable

pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'separable spatial (ISHT)');


%%%% non-separable tests

% Generate a non-separable (in theta and phi) mesh by euler rotation
[theta,phi]=ndgrid(tv,pv); % start with separable mesh

% spherical to cartesians on unit sphere (based on sph2cart.m)
ST=sin(theta);
X=ST.*cos(phi);
Y=ST.*sin(phi);
Z=cos(theta);

% rotate in cartesians
R=RZRYRZdeg(-30,80,0); % rotate matrix
Rxyz=R*[X(:)'; Y(:)'; Z(:)']; % rotate points in mesh
X(:)=Rxyz(1,:)'; % poke back into position in matrix
Y(:)=Rxyz(2,:)';
Z(:)=Rxyz(3,:)';

% cartesians to spherical on unit sphere (based on cart2sph.m)
theta=atan2(hypot(X,Y),Z);
phi=atan2(Y,X);


%%%% TEST non-separable sphereHarm

[Ylm,~,~]=sphereHarm(l,m,theta,phi);

pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-separable sphereHarm');


%%%% TEST non-separable spatial (ISHT)

[Ylm,~,~]=spatial(w,theta,phi,1); % non-separable

pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-separable spatial (ISHT)');

function doPlot(Ylm,theta,phi,pos,name)
maxY=max(abs(Ylm(:)));
Ylm=Ylm/maxY;
figure
bump_height=0.8; ref_sphere=1.0; plottype=2; % real
s=spatialPlot(Ylm,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='black'; % 'none'
fig=gcf;
fig.Position(1)=pos;
fig.Position(2)=50+pos/5;
fig.Name=name;