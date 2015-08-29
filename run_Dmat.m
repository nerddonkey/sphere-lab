function run_Dmat
 	D=Dmat_rect(7,pi/6,pi/4,pi/8,3*pi/8);

	disp('Size of D matrix')
	disp(size(D))

% 	region is whole of sphere so D should be the identity
% 	D=Dmat_rect(2,0,pi,0,2*pi);
	
	[V,E] = eigs(D);
	disp('Largest Eigenvalue')
	disp(E(1,1))

	w=V(:,1);
	
	check=bsxfun(@rdivide,D*w,w)
end
