function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     C1 and R1: the first camera pose, size (3 x 1) and (3 x 3) respectively
%     C2 and R2: the second camera pose, size (3 x 1) and (3 x 3) respectively
%     C3 and R3: the third camera pose, size (3 x 1) and (3 x 3) respectively
%     x1, x2, and x3: N × 2 matrices whose row represents correspondence
%                     between the first, second, and third images where N is 
%                     the number of correspondences.
%     X0: size (Nx3) linearly triangulated points
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations
X = []; 
for i = 1:length(x1)
  x_1 = x1(i,:);
  x_2 = x2(i,:);
  x_3 = x3(i,:);
  X_0 = X0(i,:);
  
  X_ = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x_1, x_2, x_3, X_0)
  X = [X; X_];
end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
  % concatenate correspondences
  b = [x1, x2, x3]';
  
  % projection points
  uv_1 = K*R1*(X0'-C1);
  uv_2 = K*R2*(X0'-C2);
  uv_3 = K*R3*(X0'-C3);
  uv_1 = [uv_1(1), uv_1(2)]/uv_1(3);
  uv_2 = [uv_2(1), uv_2(2)]/uv_2(3);
  uv_3 = [uv_3(1), uv_3(2)]/uv_3(3);
  fX = [uv_1, uv_2, uv_3]';
  
  % Jacobian
  J1 = Jacobian_Triangulation(C1, R1, K, X0);
  J2 = Jacobian_Triangulation(C2, R2, K, X0);
  J3 = Jacobian_Triangulation(C3, R3, K, X0);
  J = [J1', J2', J3']';
  
  % calculate step
  delta_X = inv(J'*J) * J' * (b-fX);
  
  % update X
  X = X0 + reshape(delta_X,1,3);
end

function J = Jacobian_Triangulation(C, R, K, X)
  uvw = K*R*(X'-C);
  ui = uvw(1);
  vi = uvw(2);
  wi = uvw(3);
  % partial derivatives
  du_X = [K(1,1), K(1,3)]*[[R(1,1);R(3,1)], [R(1,2);R(3,2)], [R(1,3);R(3,3)]];
  dv_X = [K(2,2), K(2,3)]*[[R(2,1);R(3,1)], [R(2,2);R(3,2)], [R(2,3);R(3,3)]];
  dw_X = R(3,:);
  
  J = [(wi*du_X - ui*dw_X)/(wi*wi); (wi*dv_X - vi*dw_X)/(wi*wi)];
end

end