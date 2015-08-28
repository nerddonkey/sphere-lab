function [f,theta,phi]=spatial_globe(L_max,ntt,npp,filename)
	Q_max=round((L_max+2)*(L_max+1)/2); % number of m>=0 spherical harmonics
	N_max=(L_max+1)^2; % total number of spherical harmonics
	disp(sprintf('\n* Numbers of Coefficients\n'))
	disp(' L_max  Q_max  N_max')
	disp(sprintf('%6d %6d %6d', L_max, Q_max, N_max))

	P_max=min(15,Q_max); % number of preview real SH file lines

	% READ IN A FILE OF REAL SH COEFFICIENTS
	sizeA=[4 Q_max];
	A=zeros(sizeA); % allocate
	fid=fopen(filename,'r');
	A=fscanf(fid,'%d %d %f %f',sizeA);
	C=A(3,:); % the cos real SH coefficient terms
	S=A(4,:); % the sin real SH coefficient terms
	disp(sprintf('\n* Top lines used from data file: %s\n', filename))
	disp(sprintf('%6d %6d %22.13e %22.13e \n',A(:,[1:P_max])))
	disp(sprintf('* Bottom lines used from data file: %s\n', filename))
	disp(sprintf('%6d %6d %22.13e %22.13e \n',A(:,[Q_max-P_max:Q_max])))
	fclose(fid);
	clear A;

	% convert the real SH from file to complex SH coefficients
	C(1)=0; % remove DC or radius term
	S(1)=0; % remove DC or radius term
 	w=zeros(N_max,1); % allocate for complex SH coefficients
	for l=[0:L_max]
		for m=[0:l]
			n=l*(l+1)+m; % (7.39) corresponding to l,m (m>=0) [output]
			n1=l*(l+1)-m; % n corresponding to l,-m (m>=0) [output]
			q=round(l*(l+1)/2+m); % real SH index [input]
			if m==0
				w(n+1)=C(q+1)-1j*S(q+1);
			else
				w(n+1)=sqrt(0.5)*(C(q+1)-1j*S(q+1));
				w(n1+1)=(-1)^m*conj(w(n+1));
			end
		end
	end

	% display some of the converted complex SH coefficients
	disp(sprintf('* Some of the converted complex SH coefficients\n'))
	R_max=min(25,L_max); % equivalent number of preview complex SHs
	disp('     n      l      m               Complex SH Coefficients');
	for n=[1:R_max]
		l=floor(sqrt(n-1)); m=n-1-l*(l+1); %(7.40)
		disp(sprintf('%6d %6d %6d %22.13e + %22.13e j', n, l, m, real(w(n)), imag(w(n))))
	end

	% evaluate the spatial function on mesh
	tt=linspace(0,pi,ntt);
	pp=linspace(0,2*pi,npp);
	[f,theta,phi]=spatial(w,tt,pp);
end
