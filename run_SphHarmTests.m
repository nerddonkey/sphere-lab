function run_SphHarmTests2
%run_SphereHarm tests various configurations of teh computation of
% spherical harmonics on separable and non-separable meshes

close all

l=7;
m=-2;

% pick delta spherical harmonic coefficients for ISHT
n=l*(l+1)+m; % (7.39)
N_tot=ceil(sqrt(n+1))^2; % max number of coefficients
w=zeros(N_tot,1); % allocate for complex SH coefficients
w(n+1)=1;

deginc=5.0;
tv=(0:deginc:180)*pi/180;
pv=(0:deginc:360)*pi/180;

pos=50;pinc=150;


%%%% TEST separable sphHarmGrid
tic
[Ylm,theta,phi]=sphHarmGrid(l,m,tv,pv);
fprintf('sphHarmGrid: %.4f seconds\n',toc);
doPlot(Ylm,theta,phi,pos,'separable sphHarmGrid');


%%%% TEST separable sphereHarmBankGrid
tic
[YB,theta,phi]=sphHarmBankGrid(l,tv,pv);
Ylm=YB{abs(m)+1};
if m<0
	Ylm=(-1)^m*conj(Ylm);
end
fprintf('sphHarmBankGrid: %.4f seconds\n',toc);
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'separable sphHarmBankmGrid');


%%%% TEST separable ishtGrid (ISHT)
tic
[Ylm,theta,phi]=ishtGrid(w,tv,pv,1);
fprintf('ishtGrid: %.4f seconds\n',toc);
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'separable ishtGrid');


%%%% TEST irregular separable sphHarmGrid
tic
r=randspace(0,70,18);r(1)=0;r(end+1)=18;tv1=r*pi/18;
r=randspace(0,140,36);r(1)=0;r(end+1)=36;pv1=r*pi/18;
[Ylm,theta,phi]=sphHarmGrid(l,m,tv1,pv1);
fprintf('sphHarmGrid: %.4f seconds\n',toc);
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'irregular separable sphHarmGrid');


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


%%%% TEST non-separable sphHarm
tic
Ylm=sphHarm(l,m,theta,phi);
fprintf('sphHarm: %.4f seconds\n',toc)
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-separable sphHarm');


%%%% TEST non-separable isht (inverse SHT)

tic
F=isht(w,theta,phi,1); % non-separable
fprintf('isht: %.4f seconds\n',toc)
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-separable isht');

function doPlot(Ylm,theta,phi,pos,name)
maxY=max(abs(Ylm(:)));
Ylm=Ylm/maxY;
figure
bump_height=0.8; ref_sphere=1.0; plottype=2; % real
s=spatial_plot(Ylm,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='black'; % 'none'
fig=gcf;
fig.Position(1)=pos;
fig.Position(2)=50+pos/5;
fig.Name=name;
