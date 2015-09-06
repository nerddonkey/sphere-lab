function run_SphereHarm
%run_SphereHarm tests sphereHarm

l=40;
m=-1;
deginc=0.25;
tt=[0:deginc:180]*pi/180;
pp=[0:deginc:360]*pi/180;

%%%% TEST sphereHarm

[Y,theta,phi]=sphereHarm(l,m,tt,pp);

% normalize entries to interval [-1.0,1.0]
maxY=max(abs(Y(:)));
Y=Y/maxY;

bump_height=0.8; ref_sphere=1.0; plottype=2; % real
s=spatial_plot(Y,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='none'; % no lines

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
s=spatial_plot(Y,theta,phi,bump_height,ref_sphere,plottype);
s.EdgeColor='none'; % no lines
