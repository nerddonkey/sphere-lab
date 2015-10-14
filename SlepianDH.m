function [D]=SlepianDH(Lmax,thetaVecR,phiVecR,maskR,useProgressBar)
%SlepianDH generate the Hermitian Slepian D matrix for mesh region
%  \int_R I_R(\unit{x}) Y_p^q(\unit{x}) \conj{Y_l^m(\unit{x})} ds(\unit{x})
%  Lmax is the highest non-zero degree
%  maskR and phiVecR define the integration region R
%  maskR irregular region mask for I_R(\unit{x}) with 0-1 matrix on maskR-phiVecR grid

addpath code/ % progressbar

%% Defaults
if nargin<5
	useProgressBar=false;
end
if nargin < 4
	maskR=false; % no mask; trapSphereMaskedR will bypass applying mask
end

%% Pre-allocate D
Ntot=(Lmax+1)^2;
D=zeros(Ntot,Ntot);

%% Number of elements to compute for Hermitian (the remaining Ntot*(Ntot-1)/2 can be inferred)
total=Ntot*(Ntot+1)/2;
if useProgressBar; count=0; progressbar('SlepianDH overall','inner loop'); end;

%% Populate the Hermitian D using (8.27)
 for r=1:Ntot % rows of D
	if useProgressBar && count~=0; progressbar([],0); end;
	l=floor(sqrt(r-1)); m=r-1-l*(l+1); % (7.40) here n=r-1
	f=sphHarmGrid(l,m,thetaVecR,phiVecR);
	D(r,r)=trapSphereMaskedR(f.*conj(f),thetaVecR,phiVecR,maskR); % diagonal element
	if useProgressBar; count=count+1; progressbar(count/total,[]); end
	for c=r+1:Ntot % columns of D
		p=floor(sqrt(c-1)); q=c-1-p*(p+1); % (7.40) here n=c-1
		g=sphHarmGrid(p,q,thetaVecR,phiVecR);
		D(r,c)=trapSphereMaskedR(g.*conj(f),thetaVecR,phiVecR,maskR); % (8.27)
		D(c,r)=conj(D(r,c)); % Hermitian symmetric element
		if useProgressBar; count=count+1; progressbar(count/total,(c-r)/(Ntot-r)); end
	end
	if useProgressBar; progressbar(count/total,[]); end;
end
