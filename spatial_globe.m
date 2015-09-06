function [F,theta,phi]=spatial_globe(L_max,ntt,npp,filename)
%spatial_globe reads files in the l m C S real spherical harmonic format
% and generates the spatial function
% For data see https://geodesy.curtin.edu.au/research/models/Earth2012/
% or http://www.ipgp.fr/~wieczor/SH/SH.html

Q_tot=round((L_max+2)*(L_max+1)/2); % number of m>=0 spherical harmonics
N_tot=(L_max+1)^2; % total number of spherical harmonics
fprintf('\n@@ Numbers of Coefficients\n\n   L_max    Q_max    N_max\n')
fprintf('%8d %8d %8d', L_max, Q_tot, N_tot)

P_max=min(15,Q_tot); % number of preview real SH file lines

% READ IN A FILE OF REAL SH COEFFICIENTS
sizeA=[4 Q_tot]; % have to read in transpose (4 rows)
A=zeros(sizeA); % allocate
fid=fopen(filename,'r');
A=fscanf(fid,'%d %d %f %f',sizeA);
C=A(3,:); % the cos real SH coefficient terms
S=A(4,:); % the sin real SH coefficient terms
fprintf('\n\n@@ Top lines used from data file: %s\n\n', filename)
fprintf('%6d %6d %22.13e %22.13e \n',A(:,[1:P_max]))
fprintf('\n@@ Bottom lines used from data file: %s\n\n', filename)
fprintf('%6d %6d %22.13e %22.13e \n',A(:,[Q_tot-P_max:Q_tot]))
fclose(fid);
clear A;

% convert the real SH from file to complex SH coefficients
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
			w(n1+1)=(-1)^m*conj(w(n+1)); % Y_l^{-}m coefficient (m>0)
		end
	end
end

% display some of the converted complex SH coefficients
fprintf('\n@@ Some of the converted complex SH coefficients\n\n')
R_max=min(25,L_max); % equivalent number of preview complex SHs
disp('     n      l      m               Complex SH Coefficients');
for n=[1:R_max]
	l=floor(sqrt(n-1)); m=n-1-l*(l+1); %(7.40)
	fprintf('%6d %6d %6d %21.13e%+21.13ej\n',n,l,m,real(w(n)),imag(w(n)))
end

% evaluate the spatial function on mesh
if 1 % Australia
	tt=linspace(90,140,ntt)*pi/180;
	pp=linspace(285,345,npp)*pi/180;
else % whole world
	tt=linspace(0,180,ntt)*pi/180;
	pp=linspace(0,360,npp)*pi/180;
end
[F,theta,phi]=spatial(w,tt,pp,1);
maxF=max(abs(F(:)));
F=F/maxF; % normalize entries to interval [-1.0,1.0]

