function [ H ] = est_homography(video_pts, logo_pts)
% est_homography estimates the homography to transform each of the
% video_pts into the logo_pts
% Inputs:
%     video_pts: a 4x2 matrix of corner points in the video
%     logo_pts: a 4x2 matrix of logo points that correspond to video_pts
% Outputs:
%     H: a 3x3 homography matrix such that logo_pts ~ H*video_pts
% Written for the University of Pennsylvania's Robotics:Perception course

% YOUR CODE HERE
A = [];
for i = 1:4
  x = [video_pts(i,1), video_pts(i,2)];
  x_p = [logo_pts(i,1), logo_pts(i,2)];

  a_x = [-x(1), -x(2), -1, 0, 0, 0, x(1)*x_p(1), x(2)*x_p(1), x_p(1)];
  a_y = [0, 0, 0, -x(1), -x(2), -1, x(1)*x_p(2), x(2)*x_p(2), x_p(2)];
  A = vertcat(A, a_x, a_y);
end
[U, S, V] = svd(A);
h = V(:,end);
H = reshape(h, [3,3])';

end

