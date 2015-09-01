function [D]= SlepianD(L_max,tt,pp)
%SlepianD generate the Slepian D matrix for rect region
% to do: use Hermitian property to roughly halve the computation
	N_max=(L_max+1)^2;
	D=zeros(N_max,N_max); % pre-allocate D
	progressbar('row index','column index')
	for r=1:N_max % rows
		progressbar([],0)
		for c=1:N_max % columns
			D(r,c)=SlepianDrc(r,c,tt,pp);
			progressbar((r-1)/N_max+c/N_max^2,c/N_max)
		end
		progressbar(r/N_max,[])
	end
	D=(D+D')/2; % make perfect Hermitian
end
