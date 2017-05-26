function run_Rotsymtest
%run_Rotsymtest render a superposition of symmetric functions

close all
clear ll w_eta w_mu w

%% range of degree
L=30; % maximum (inclusive) degree
N_tot=(L+1)^2; % max number of SH coefficients
ll=0:L; % index 1,2,...,L+1 convenience vector


%% 1: eta - north pole component

w_eta=zeros(N_tot,1); % 1D coeffs initially filled with zeros

kappa_eta=30; % pointiness
lambda_eta=exp(-ll.*(ll+1)/(2*kappa_eta)) % Gauss-Weierstrass kernel

% use north pole formula to determine the SH coefficients of eta component
for l=0:L
	m=0; % just a reminder
	n=l*(l+1)+m; % (7.39)
	w_eta(n+1)=sqrt((2*l+1)/(4*pi)) * lambda_eta(l+1);
end


%% 2: mu - arbitrary angle component

w_mu=zeros(N_tot,1); % 1D coeffs initially filled with zeros

% mean direction (anything)
theta_mu=0.55*pi;
phi_mu=3.5*pi;

kappa_mu=140; % pointiness
lambda_mu=exp(-ll.*(ll+1)/(2*kappa_mu)) % Gauss-Weierstrass kernel

% use formula to determine the SH coefficients of mu component
for l=0:L
	for m=-l:l
		n=l*(l+1)+m; % (7.39)
		w_mu(n+1)=lambda_mu(l+1)*conj(sphHarm(l,m,theta_mu,phi_mu)); % scalar angles here
	end
end

%% SH coefficients of superposition
w=2.0*w_eta+1.0*w_mu; % some superposition

%% determine spatial samples on mesh
deginc=2.0;
tv=(0:deginc:180)*pi/180;
pv=(-30:deginc:330)*pi/180;	
[w_spatial,theta,phi]=ishtRectGrid(w,tv,pv,true);

%% do plot
doThePlot(w_spatial,theta,phi,'rotationally symmetric function plot');


function doThePlot(Ylm,theta,phi,name)
%doThePlot - convenience function to plot the spatial form
maxY=max(abs(Ylm(:)));
Ylm=Ylm/maxY; % normalize
figure
bump_height=0.8; ref_sphere=1.0;
plottype=2; % real
s=spatialPlot(Ylm,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='black'; % 'none'
fig=gcf;
fig.Name=name;
fig.Position(3)=600;
fig.Position(4)=500;
alpha 0.8
