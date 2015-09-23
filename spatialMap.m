function F=spatialMap(G,plottype,doScale)
%spatialMap

%% default arguments
if nargin<2
	plottype=1;
end
if nargin<3
	doScale=0;
end

%% perform mapping
switch plottype
	case 0
		F=abs(G).^2;
	case 1
		F=abs(G);
	case 2
		F=real(G);
	case 3
		F=imag(G);
	otherwise
		F=real(G);
end

%% scale entries to interval [-1.0,1.0]
if doScale~=0
 	maxF=max(abs(F(:)));
 	F=F/maxF;
end
