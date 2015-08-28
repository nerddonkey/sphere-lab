function [D]=run_Dmat
% 	D=Dmat_rect(1,pi/6,pi/4,pi/8,3*pi/8);
	D=Dmat_rect(1,0,pi,0,2*pi);
	eig(D)
end
