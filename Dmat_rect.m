function [D]= Dmat_rect(L,t1,t2,p1,p2)
%Dmat generate the Slepian D matrix for rect region
	N_max=(L+1)^2;
	D=zeros(N_max,N_max); % pre-allocate D

	delta=0.01; % roughly 0.5730 degrees increment or less
	ntt=ceil((t2-t1)/delta)+1;
	npp=ceil((p2-p1)/delta)+1;

	tt=linspace(t1,t2,ntt);
	pp=linspace(p1,p2,npp);

	for r=1:N_max % rows
		for c=1:N_max % columns
			D(r,c)=herm(r,c,tt,pp);
		end
	end
end
