function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
%% LinearTriangulation
% Find 3D positions of the point correspondences using the relative
% position of one camera from another
% Inputs:
%     C1 - size (3 x 1) first camera center in the world coordinate
%     R1 - size (3 x 3) rotation of the first camera pose
%     C2 - size (3 x 1) second camera center in the world coordinate
%     R2 - size (3 x 3) rotation of the second camera pose
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Outputs: 
%     X - size (N x 3) matrix whos rows represent the 3D triangulated
%       points

% Camera matrices
P1 = K*[R1, -R1*C1]; % t1 = -R1*C1 is position of the world coordinate w.r.t camera center
P2 = K*[R2, -R2*C2]; % C2 is position of cam2 in the world coordinate w.r.t camera center
X = [];
% Correspondences
for i = 1:length(x1)
  x_1 = [x1(i,:)';1];
  x_2 = [x2(i,:)';1];
  skew1 = Vec2Skew(x_1);
  skew2 = Vec2Skew(x_2);
  A = [skew1 * P1; skew2 * P2];
  % Solve
  [~,~,v] = svd(A);
  X_ = v(:,end)/v(end,end);
  X = [X; X_(1:3)'];
end

