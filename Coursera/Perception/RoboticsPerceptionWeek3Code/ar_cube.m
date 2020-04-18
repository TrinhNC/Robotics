function [proj_points, t, R] = ar_cube(H,render_points,K)
%% ar_cube
% Estimate your position and orientation with respect to a set of 4 points on the ground
% Inputs:
%    H - the computed homography from the corners in the image
%    render_points - size (N x 3) matrix of world points to project
%    K - size (3 x 3) calibration matrix for the camera
% Outputs: 
%    proj_points - size (N x 2) matrix of the projected points in pixel
%      coordinates
%    t - size (3 x 1) vector of the translation of the transformation
%    R - size (3 x 3) matrix of the rotation of the transformation
% Written by Stephen Phillips for the Coursera Robotics:Perception course

% YOUR CODE HERE: Extract the pose from the homography
if H(9) < 0
    % enforce the z of t to be positive
    H = H * [1,0,0;0,1,0;0,0,-1];
end
h1 = H(:,1);
h2 = H(:,2);
h3 = H(:,3);
hx = cross(h1, h2);
R_prime = [h1, h2, hx];
[U,S,V] = svd(R_prime);
R = U * [[1, 0, 0];[0, 1, 0];[0, 0, det(U*V)]] * V;
t = h3/norm(h1);
% YOUR CODE HERE: Project the points using the pose
Xc = K * (R*render_points' + t);
zc = abs(Xc(3,:));
proj_points = Xc./zc;
proj_points = proj_points(1:2, :)';
end
