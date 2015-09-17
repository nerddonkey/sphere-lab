function run_SphereHarmTests
%run_SphereHarm tests sphereHarm

close all

l=5;
m=-1;
deginc=5.0;
tt=(0:deginc:180)*pi/180;
pp=(0:deginc:360)*pi/180;

%%%% TEST sphereHarm

[Y,theta,phi]=sphereHarm(l,m,tt,pp);

% normalize entries to interval [-1.0,1.0]
maxY=max(abs(Y(:)));
Y=Y/maxY;

bump_height=0.8; ref_sphere=1.0; plottype=2; % real
s=spatial_plot(Y,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='black'; % 'none'
fig=gcf;
fig.Position(1)=0;
fig.Position(2)=100;

figure

%%%% TEST sphereHarmBank

[YlBank,theta,phi]=sphereHarmBank(l,tt,pp);

Ylm=YlBank{abs(m)+1};
if m<0
	Ylm=(-1)^m*conj(Ylm);
end

maxY=max(abs(Ylm(:)));
Ylm=Ylm/maxY;
bump_height=0.8; ref_sphere=1.0; plottype=2; % real
s=spatial_plot(Ylm,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='black'; %'none'
fig=gcf;
fig.Position(1)=400;
fig.Position(2)=100;

figure

%%%% non-separable test

% Generate a non-separable (in theta and phi) mesh by euler rotation
[theta,phi]=ndgrid(tt,pp); % start with separable mesh

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

% pick delta spherical harmonic coefficients
n=l*(l+1)+m; % (7.39)
N_tot=ceil(sqrt(n+1))^2; % max number of coefficients
w=zeros(N_tot,1); % allocate for complex SH coefficients
w(n+1)=1;

[Ylm,theta,phi]=spatial(w,theta,phi,1); % non-separable

maxF=max(abs(Ylm(:)));
Ylm=Ylm/maxF;
bump_height=0.8; ref_sphere=1.0; plottype=2; % real
s=spatial_plot(Ylm,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='black'; % 'none'
fig=gcf;
fig.Position(1)=800;
fig.Position(2)=100;

