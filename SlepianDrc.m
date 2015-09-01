function [y,theta,phi]=SlepianDrc(r,c,tt,pp)
%SlepianDrc compute D_{r,c} in (8.27) where r={l,m} and c={p,q}
% here r=row and c=column
% tt and pp should be uniform vectors (at this stage)
% To do: non-uniform grid case; use trapz
	[f,theta,phi]=spatial_shn(r,tt,pp); % Y_r or Y_l^m in (8.27)
	[g,~,~]=spatial_shn(c,tt,pp); % Y_c or Y_p^q in (8.27)
	f=g.*conj(f); % Y_p^q \conj{Y_l^m} in (8.27)
	if 0
		f=bsxfun(@times,f,sin(tt(:))); % Y_p^q \conj{Y_l^m} sin(theta) in (8.27)
		dttdpp=(tt(2)-tt(1))*(pp(2)-pp(1)); % dtheta dphi
		y=sum(f(:))*dttdpp; % sum to approximate integration (8.27)
	else
		y=trapSphereR(f,tt,pp);
	end
end
