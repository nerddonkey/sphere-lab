function [tv,pv,R_mask,R_theta,R_phi]=ausRegion(deginc,plotme)
%ausRegion
% grabs MATLAB's built-in world coastline data: coast
% pulls out mainland Australia and Tasmania
% show how to combine to a single "polygon" representation: R_theta and R_phi
% creates a separable mesh [theta,phi] containing the polygon
% creates R_mask on the (separable) mesh which is 1 inside the polygon and 0 outside
% optionally plots the results to verify it makes sense

if nargin<2
	plotme=0;
end
if nargin<1
	deginc=1.0; % default stepsize in degrees
end

load coast % returns column vectors long and lat in degrees
R_main=8296:8604; % mainland Australia; 8604 is the same as 8296
R_tas=8623:8645; % Tasmania; 8645 is the same as 8623

% define the region R (in degrees) as a union
R_theta=[(90-lat(R_main))' NaN (90-lat(R_tas))'];
R_phi=[long(R_main)' NaN long(R_tas)'];

% find whole-degree enlarged bounding rectangle on region
min_theta=floor(min(R_theta));
max_theta=ceil(max(R_theta));
min_phi=floor(min(R_phi));
max_phi=ceil(max(R_phi));

% determine covering mesh with step <= deginc
numtt=ceil((max_theta-min_theta)/deginc);
numpp=ceil((max_phi-min_phi)/deginc);
tv=linspace(min_theta,max_theta,numtt)*pi/180;
pv=linspace(min_phi,max_phi,numpp)*pi/180;

[theta,phi]=ndgrid(tv,pv);

% change the coastline data from degrees to radians
R_theta=R_theta*pi/180;
R_phi=R_phi*pi/180;

% 0-1 mask for whole-degree enlarged bounding rectangle
R_mask=inpolygon(phi,theta,R_phi,R_theta);
% usage: I=trapSphereMaskedR(f,tt,pp,R_mask) for some f

if plotme~=0 % plot on sphere to see if makes sense
	hold on
	% coastline
	aus_coast=[sin(R_theta).*cos(R_phi); sin(R_theta).*sin(R_phi); cos(R_theta)];
	idx = any([~isnan(R_theta);~isnan(R_theta)],1); %This is to not plot if either x or y has a nan
	fill3(aus_coast(1,idx),aus_coast(2,idx),aus_coast(3,idx),'yellow','EdgeColor','None');
	% points inside
	pQ=phi(R_mask);
	tQ=theta(R_mask);
	points=[sin(tQ(:)).*cos(pQ(:)), sin(tQ(:)).*sin(pQ(:)), cos(tQ(:))]';
	plot3(points(1,:),points(2,:),points(3,:),'b.') % points inside
	% points outside
	pQ=phi(~R_mask);
	tQ=theta(~R_mask);
	points=[sin(tQ(:)).*cos(pQ(:)), sin(tQ(:)).*sin(pQ(:)), cos(tQ(:))]';
	plot3(points(1,:),points(2,:),points(3,:),'r.') % points outside

	view(-130,-30) % set viewpoint (suits Australia)
	axis off
	set(gca,'CameraViewAngle',7) % zoom in to frame Australia
	shg
	hold off
end
