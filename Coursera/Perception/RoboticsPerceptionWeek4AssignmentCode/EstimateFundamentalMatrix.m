function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

% Construct matrix A
A = [];
for i = 1:length(x1)
  u1 = x1(i,1);
  v1 = x1(i,2);
  u2 = x2(i,1);
  v2 = x2(i,2);
  a = [u1*u2, u1*v2, u1, v1*u2, v1*v2, v1, u2, v2, 1];
  A = [A;a];
end

% Solving linear homogeneous equations via SVD
[u,s,v] = svd(A);
F = v(:,end);
%F = v(:,end)/v(end,end);
F = reshape(F,3,3)';

% Apply rank constraint (rank(F) must be 2)
[u,s,v] = svd(F);
s(end,end) = 0;
F = u*s*v';
F = F/norm(F);
end
