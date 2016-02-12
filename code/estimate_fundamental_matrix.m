% Fundamental Matrix Stencil Code
% CS 4495 / 6476: Computer Vision, Georgia Tech
% Written by Henry Hu

% Returns the camera center matrix for a given projection matrix

% 'Points_a' is nx2 matrix of 2D coordinate of points on Image A
% 'Points_b' is nx2 matrix of 2D coordinate of points on Image B
% 'F_matrix' is 3x3 fundamental matrix

% Try to implement this function as efficiently as possible. It will be
% called repeatly for part III of the project

function [ F_matrix ] = estimate_fundamental_matrix(Points_a,Points_b)

%%%%%%%%%%%%%%%%
% Your code here
ua=Points_a(:,1);
va=Points_a(:,2);
ub=Points_b(:,1);
vb=Points_b(:,2);

L=size(Points_a,1);
% 8point
l1=ones(L,1);
%A=[ua.*ub, ua.*vb, ua, va.*ub, va.*vb, va, ub, vb, l1];
A=[ub.*ua, ub.*va, ub, vb.*ua, vb.*va, vb, ua, va, l1];
[U,S,V]=svd(A);
f=V(:,end);
F=reshape(f,[3 3])';

[UU,SS,VV]=svd(F);
SS(3,3)=0;
F_matrix=UU*SS*VV';

%%%%%%%%%%%%%%%%

%This is an intentionally incorrect Fundamental matrix placeholder
%F_matrix = [0  0     -.0004; ...
%            0  0      .0032; ...
%            0 -0.0044 .1034];
        
end

