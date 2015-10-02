function [F,theta,phi,L_max,radius]=ishtFromShapeFile(L_req,ntv,npv,filename,doScale)
%ishtFromShapeFile inverse spherical harmonic transform (spatial function)
%  from *.shape real spherical harmonic coefficient planetary data
%
%  Reads *.shape files in the l m C S real spherical harmonic coefficient
%  format, converts to complex spherical harmonic coefficients and
%  generates the spatial function on the sphere.
%
%  Generates a whole sphere rectangular mesh (theta,phi) with ntv samples
%  in theta and npv samples in phi.  The spatial function (F) is evaluated
%  on this mesh.
%
%  For sample shape planetary data see:
%    https://geodesy.curtin.edu.au/research/models/Earth2012/
%    http://www.ipgp.fr/~wieczor/SH/SH.html

%% Attempt to read *.shape file
[fid,errMessage]=fopen(filename,'r');
if fid==-1
	error(['@@ ishtFromShapeFile Error: "' filename '" - ' errMessage]);
end

if nargin<5
	doScale=1;
end

%% Try to read for requested degree
Q_req=round((L_req+2)*(L_req+1)/2); % number of m>=0 spherical harmonics

%% Read in real SH coefficients from shape file
sizeA=[4,Q_req]; % have to read in transpose (4 rows)
A=fscanf(fid,'%d %d %f %f',sizeA);
C=A(3,:); % the cos real SH coefficient terms
S=A(4,:); % the sin real SH coefficient terms
radius=C(1); % usually but not always the mean radius
L_max=A(1,end); % this will be the min of L_req and the actual degrees available

Q_tot=round((L_max+2)*(L_max+1)/2); % number of m>=0 spherical harmonics
N_tot=(L_max+1)^2; % total number of spherical harmonics
fprintf('\n@@ Numbers of Coefficients\n\n   L_max    Q_max    N_max\n')
fprintf('%8d %8d %8d', L_max, Q_tot, N_tot)

%% Display top and bottom of data to be used
P_max=min(15,Q_tot); % number of preview real SH file lines
fprintf('\n\n@@ Top lines used from data file: %s\n\n', filename)
fprintf('%6d %6d %22.13e %22.13e \n',A(:,[1:P_max]))
fprintf('\n@@ Bottom lines used from data file: %s\n\n', filename)
fprintf('%6d %6d %22.13e %22.13e \n',A(:,[Q_tot-P_max:Q_tot]))
fclose(fid);
clear A;

%% Convert the real SH coefficients to complex SH coefficients
C(1)=0; % remove DC or radius term
S(1)=0; % remove DC or radius term (should be zero anyway)
w=zeros(N_tot,1); % allocate for complex SH coefficients
for l=[0:L_max]
	for m=[0:l]
		n=l*(l+1)+m; % (7.39) corresponding to l,m (m>=0) [output]
		n1=l*(l+1)-m; % n corresponding to l,-m (m>=0) [output]
		q=round(l*(l+1)/2+m); % real SH index [input]
		if m==0
			w(n+1)=C(q+1)-1j*S(q+1); % Y_l^0 coefficient (m=0)
		else
			w(n+1)=sqrt(0.5)*(C(q+1)-1j*S(q+1)); % Y_l^m coefficient (m>0)
			w(n1+1)=(-1)^m*conj(w(n+1)); % Y_l^{-m} coefficient (m>0)
		end
	end
end

%% Display some of the converted complex SH coefficients
fprintf('\n@@ Some of the converted complex SH coefficients\n\n')
fprintf('     n      l      m               Complex SH Coefficients\n');
R_max=min(25,L_max); % equivalent number of preview complex SHs
for n=[0:R_max-1]
	l=floor(sqrt(n)); m=n-l*(l+1); % (7.40)
	fprintf('%6d %6d %6d %21.13e%+21.13ej\n',n,l,m,real(w(n+1)),imag(w(n+1)))
end

%% Define separable mesh vectors
if 0 % Australia
	ntv=ntv*5/18;
	tv=linspace(90,140,ntv)*pi/180;
	npv=npv/6;
	pv=linspace(285,345,npv)*pi/180;
else % whole world
	tv=linspace(0,180,ntv)*pi/180;
	pv=linspace(0,360,npv)*pi/180;
end

%% Perform the ISHT from w to F on separable mesh [theta,phi]=ndgrid(tv,pv)
% slower [F,theta,phi]=ishtGrid2(w,tv,pv,1);
[F,theta,phi]=ishtRectGrid(w,tv,pv,true,true); % useProgreebar and doReal

%% Normalize entries to interval [-1.0,1.0]
if doScale
	maxF=max(abs(F(:)));
	F=F/maxF;
end
