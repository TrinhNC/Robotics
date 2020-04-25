function E = EssentialMatrixFromFundamentalMatrix(F,K)
%% EssentialMatrixFromFundamentalMatrix
% Use the camera calibration matrix to esimate the Essential matrix
% Inputs:
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     F - size (3 x 3) fundamental matrix from EstimateFundamentalMatrix
% Outputs:
%     E - size (3 x 3) Essential matrix with singular values (1,1,0)

E = K' * F * K;
[u,s,v] = svd(E);
s = [1,0,0;0,1,0;0,0,0];
E = u * s * v';
E = E/norm(E); % need to normalize because we are in homogeneous coordinate