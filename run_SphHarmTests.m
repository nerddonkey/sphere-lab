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

deginc=2.0;
tv=(0:deginc:180)*pi/180;
pv=(0:deginc:360)*pi/180;

pos=50;pinc=150;


%%%% TEST separable sphHarmGrid
tic
[Ylm,theta,phi]=sphHarmGrid(l,m,tv,pv);
fprintf('uniform separable sphHarmGrid: %.4f seconds\n',toc);
doPlot(Ylm,theta,phi,pos,'uniform separable sphHarmGrid');


% %%%% TEST separable sphereHarmBankGrid
% tic
% [YB,theta,phi]=sphHarmBankGrid(l,tv,pv);
% Ylm=YB{abs(m)+1};
% if m<0
% 	Ylm=(-1)^m*conj(Ylm);
% end
% fprintf('sphHarmBankGrid: %.4f seconds\n',toc);
% pos=pos+pinc;
% doPlot(Ylm,theta,phi,pos,'separable sphHarmBankmGrid');


%%%% TEST separable ishtRectGrid (ISHT)
tic
[Ylm,theta,phi]=ishtRectGrid(w,tv,pv,1);
fprintf('uniform separable ishtGrid2: %.4f seconds\n',toc);
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'uniform separable ishtRectGrid');


%%%% TEST non-uniform separable sphHarmGrid
tic
r=randspace(0,70,18);r(1)=0;r(end+1)=18;tv1=r*pi/18;
r=randspace(0,140,36);r(1)=0;r(end+1)=36;pv1=r*pi/18;
[Ylm,theta,phi]=sphHarmGrid(l,m,tv1,pv1);
fprintf('non-uniform separable sphHarmGrid: %.4f seconds\n',toc);
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-uniform separable sphHarmGrid');


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
fprintf('non-separable sphHarm: %.4f seconds\n',toc)
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-separable sphHarm');


%%%% TEST non-separable isht (inverse SHT)

tic
F=isht(w,theta,phi,1); % non-separable
fprintf('non-separable isht: %.4f seconds\n',toc)
pos=pos+pinc;
doPlot(Ylm,theta,phi,pos,'non-separable isht');

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

function y = randspace(p1, n, p2, step_range, seed)
% RANDSPACE Generates a monotonically increasing sequence of randomly
%           spaced values.
%
%   y = randspace(P1, N) generates N points greater than or equal to P1
%   with intervals between the points equal to random values in [0,1].
%
%   y = randspace(P1, N, P2) cuts off the generated sequence at P2.
%   NOTE: If P2 is set, the length of the generated sequence may be less
%   than N.
%
%   y = randspace(P1, N, P2, STEP_RANGE) generates the intervals between
%   points according to STEP_RANGE.  STEP_RANGE must contain two elements
%   indicating the minimum and maximum desired intervals. The first element
%   must be less than the second.
%
%   y = randspace(P1, N, P2, STEP_RANGE, SEED) uses SEED as the state for
%   the random number generator.  SEED must be either a scalar or a 35
%   element vector. See the doc for randn for details.
%
%   Examples:
%       %generate 10 randomly spaced values starting at 0
%       y = randspace(0,10)
%
%       %generate at most 10 randomly spaced values between 1 and 6
%       y = randspace(1, 10, 6)
%
%       %generate 10 values starting at 0 with intervals between 1 and 5
%       y = randspace(0,10,[],[1,5])
%
%       %generate 10 randomly spaced values starting at 0 using one state
%       y = randspace(0,10,[],[],1)
%
%   Written by Dmitry Savransky, 18 June 2008

%check user inputs
if nargin < 2 || isempty(p1) || isempty(n)
    error('Starting point and number of points are required inputs.')
end


%if user supplied a state for rand, use it, otherwise set it to the current
%clock time.
if ~exist('seed','var') || isempty(seed)
    rand('state',sum(100*clock))
else
    if numel(seed) ~= 1 && numel(seed) ~= 35
        error('Seed must be scalar or 35-by-1.')
    end
    rand('state',seed(:));
end

%if user supplied a step range, use it, otherwise generate steps in [0 1]
if ~exist('step_range','var') || isempty(step_range)
    y = rand(1,n);
else
    if numel(step_range) ~= 2 || step_range(1) >= step_range(2)
        error('Step range must have two elements [a,b] such that a < b.')

    end
    y = step_range(1) + (step_range(2)-step_range(1)) * rand(1,n);
end

%generate sequence
y = y*triu(ones(n));
y = y+p1;

%if user supplied maximum value, use it
if exist('p2','var') && ~isempty(p2)
    y = y(y<=p2);
end
