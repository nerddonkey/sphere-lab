function run_SphHarmTests
%run_SphereHarm tests various configurations of teh computation of
% spherical harmonics on separable and non-separable meshes

close all

l=7;
m=-2;

%% Pick delta spherical harmonic coefficients for ISHT
n=l*(l+1)+m; % (7.39)
N_tot=ceil(sqrt(n+1))^2; % max number of coefficients
w=zeros(N_tot,1); % allocate for complex SH coefficients
w(n+1)=1;

deginc=2.0;
tv=(0:deginc:180)*pi/180;
pv=(-30:deginc:330)*pi/180;

pos=50;pinc=150;


%% TEST separable sphHarmGrid
tic
[Ylm,theta,phi]=sphHarmGrid(l,m,tv,pv);
fprintf('uniform separable sphHarmGrid: %.4f seconds\n',toc);
doPlot(Ylm,theta,phi,pos,'uniform separable sphHarmGrid');


%% TEST separable ishtRectGrid (ISHT)
tic
[Ylm,theta,phi]=ishtRectGrid(w,tv,pv,true);
fprintf('uniform separable ishtGrid2: %.4f seconds\n',toc);
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'uniform separable ishtRectGrid');


%% TEST non-uniform separable sphHarmGrid
tic
tv1=ranspace(0,180,90)*pi/180;
pv1=ranspace(0,360,180)*pi/180;
[Ylm,theta,phi]=sphHarmGrid(l,m,[0 tv1 pi],[0 pv1 2*pi]);
fprintf('non-uniform separable sphHarmGrid: %.4f seconds\n',toc);
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-uniform separable sphHarmGrid');


%% Test non-separable tests

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

%% Cartesians to spherical on unit sphere (based on cart2sph.m)
theta=atan2(hypot(X,Y),Z);
phi=atan2(Y,X);


%% TEST non-separable sphHarm
tic
Ylm=sphHarm(l,m,theta,phi);
fprintf('non-separable sphHarm: %.4f seconds\n',toc)
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-separable sphHarm');


%% TEST non-separable isht (inverse SHT)

tic
Ylm=isht(w,theta,phi,1); % non-separable
fprintf('non-separable isht: %.4f seconds\n',toc)
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-separable isht');


%% doPlot - convenience function to plot the spherical harmonics
function doPlot(Ylm,theta,phi,pos,name)
%doPlot
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

%% ranspace - a random version of linspace
function y=ranspace(X1,X2,N)
%ranspace RAK 28Sep15 (scalar inputs)
rr=cumsum(rand(1,N+2));
g=(X2-X1)/(rr(end)-rr(1));
y=g*(rr-rr(1))+X1;
y=y(2:end-1);
