function [thetaVecR,phiVecR,maskR,thetaR,phiR]=coastRegion(indexRange,thetaDeg,phiDeg)
%coastRegion return the coastline and mask data for region on sphere (all angles in degrees)
%
% Inputs
%   indexRange: index range in the coastline data set for the region of interest NaN delimited
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

%% Extract "island" subvectors and convert to colatitude and longitude (in degrees)

nanInd=isnan(indexRange);
indexRange(nanInd)=1; % temp dummy for NaN
thetaR=[(90-lat(indexRange))']; % row
phiR=long(indexRange)'; % row

%[thetaR,phiR]=rotatedeg(thetaR,phiR,0,155,0); % for testing

thetaR(nanInd)=NaN; % restore NaN
phiR(nanInd)=NaN; % restore NaN

%% Find whole-degree enlarged bounding rectangle on region
thetaMin=max(0.0,floor(min(thetaR)));
thetaMax=min(180.0,ceil(max(thetaR)));
phiMin=floor(min(phiR))
phiMax=ceil(max(phiR))

%% Determine covering mesh [theta,phi] with step <= thetaDeg,phiDeg
numtt=ceil((thetaMax-thetaMin)/thetaDeg);
numpp=ceil((phiMax-phiMin)/phiDeg);
thetaVecR=linspace(thetaMin,thetaMax,numtt);
phiVecR=linspace(phiMin,phiMax,numpp);
[theta,phi]=ndgrid(thetaVecR,phiVecR);

%% 0-1 mask for whole-degree enlarged bounding rectangle
maskR=inpolygon(phi,theta,phiR,thetaR);

%% Flash up the region
close
spy(maskR);
axis equal
shg
pause(2)
close

function [thetaOut,phiOut]=rotatedeg(thetaIn,phiIn,aa,bb,gg)
%% Apply SO(3)
dtr=pi/180.0;
ST=sin(thetaIn*dtr);
X=ST.*cos(phiIn*dtr);
Y=ST.*sin(phiIn*dtr);
Z=cos(thetaIn*dtr);
R=RZRYRZdeg(aa,bb,gg); % rotate matrix
Rxyz=R*[X;Y;Z]; % rotate points
X=Rxyz(1,:);
Y=Rxyz(2,:);
Z=Rxyz(3,:);
thetaOut=atan2(hypot(X,Y),Z)/dtr;
phiOut=atan2(Y,X)/dtr;

%phiR=mod(long(indexRange)',360);
