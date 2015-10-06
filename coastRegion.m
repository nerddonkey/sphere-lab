function [thetaVec,phiVec,maskR,thetaR,phiR]=coastRegion(indexRange,thetaDeg,phiDeg)
%coastRegion return the coastline and mask data for region on sphere (all angles in degrees)
%
% Inputs
%   indexRange: index range in the coastline data set for the region of interest
%   thetaDeg,phiDeg: mask and integration grid step sizes
%
% Outputs
%   thetaVec, phiVec: range vectors such that ndgrid(thetaVec,phiVec) is the close fit rectangle on region
%   maskR: is 0-1 indicator region matrix of size of outputs of ndgrid(thetaVec,phiVec)
%   thetaR,phiR: are the vectors of coastline nodes (closed curve)

if nargin<3
	phiDeg=1.0;
end
if nargin<2
	thetaDeg=1.0;
end
if nargin<1
	indexRange=8296:8604; % default to mainland Australia
end

%% Grab MATLAB's built-in world coastline data: coast
load coast % returns column vectors long and lat in degrees

%% Extract "island" subvector and convert to colatitude and longitude (in degrees)
thetaR=[(90-lat(indexRange))'];
phiR=[long(indexRange)'];

if thetaR(1)~=thetaR(end) || phiR(1)~=phiR(end)
	error('Check region index range, not explicitly a closed curve');
end

%% Find whole-degree enlarged bounding rectangle on region
thetaMin=max(0.0,floor(min(thetaR)));
thetaMax=min(180.0,ceil(max(thetaR)));
phiMin=floor(min(phiR));
phiMax=ceil(max(phiR));

%% Determine covering mesh [theta,phi] with step <= thetaDeg,phiDeg
numtt=ceil((thetaMax-thetaMin)/thetaDeg);
numpp=ceil((phiMax-phiMin)/phiDeg);
thetaVec=linspace(thetaMin,thetaMax,numtt);
phiVec=linspace(phiMin,phiMax,numpp);
[theta,phi]=ndgrid(thetaVec,phiVec);

%% 0-1 mask for whole-degree enlarged bounding rectangle
maskR=inpolygon(phi,theta,phiR,thetaR);
