function run_Orthonormal
%run_Orthonormal tests the expression for spherical harmonics
% for orthonormality.  Uses the Slepian mstrix with the region
% equal to the whole sphere.

L_max=4;
tic
deginc=0.25; % step for numerical integration
tt=(0:deginc:180)*pi/180;
pp=(0:deginc:360)*pi/180;

% interface: function [D]= SlepianDH(L_max,tt,pp,R_mask)
mask=ones(length(tt),length(pp)); % that is, no mask; just for testing
D=SlepianDH(L_max,tt,pp,mask); % uses sphereHarm; D should be identity

MSE=norm(D-eye(size(D)),'fro')/numel(D); % MSE from identity
fprintf('Size of D matrix: %d x %d\n', size(D))
fprintf('Orthonormality MSE: %8.6f\n', MSE)

D=SlepianDH(L_max,tt,pp); % uses sphereHarm; D should be identity

MSE=norm(D-eye(size(D)),'fro')/numel(D); % MSE from identity
fprintf('Size of D matrix: %d x %d\n', size(D))
fprintf('Orthonormality MSE: %8.6f\n', MSE)
toc
