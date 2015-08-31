function run_Orthonormal
%run_Orthonormal tests the expression for spherical harmonics
% for orthonormality.  Uses the Slepian mstrix with the region
% equal to the whole sphere.
   L_max=4;
   deginc=0.25; % step for numerical integration
   tt=[0:deginc:180]*pi/180;
   pp=[0:deginc:360]*pi/180;

 	D=SlepianD(L_max,tt,pp); % should be identity

	MSE=norm(D-eye(size(D)),'fro')/numel(D); % MSE from identity
	fprintf('Size of D matrix: %d x %d\n', size(D))
	fprintf('Orthonormality MSE: %8.6f\n', MSE)
end
