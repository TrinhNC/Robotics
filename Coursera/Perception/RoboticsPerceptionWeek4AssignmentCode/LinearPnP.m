function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

A = [];
for i=1:length(x)
  xc = inv(K)*[x(i,:),1]';
  xi = xc(1);
  yi = xc(2);
  Xi = X(i,:);
  A = [A;[Xi, zeros(1,3), -xi*Xi, 1, 0, -xi; zeros(1,3), Xi,  -yi*Xi, 0, 1, -yi]];
end

[u,s,v] = svd(A);
P = v(:,end);
R = reshape(P(1:9),3,3)';
t = P(10:end);

[u,d,v] = svd(R);

% clean up rotation and translation
if det(u*v')>0
  Rc = u*v';
  tc = t/d(1,1);
else
  Rc = -u*v';
  tc = -t/d(1,1);
end

R = Rc;
C = -Rc'*tc;









