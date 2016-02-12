% RANSAC Stencil Code
% CS 4495 / 6476: Computer Vision, Georgia Tech
% Written by Henry Hu

% Find the best fundamental matrix using RANSAC on potentially matching
% points

% 'matches_a' and 'matches_b' are the Nx2 coordinates of the possibly
% matching points from pic_a and pic_b. Each row is a correspondence (e.g.
% row 42 of matches_a is a point that corresponds to row 42 of matches_b.

% 'Best_Fmatrix' is the 3x3 fundamental matrix
% 'inliers_a' and 'inliers_b' are the Mx2 corresponding points (some subset
% of 'matches_a' and 'matches_b') that are inliers with respect to
% Best_Fmatrix.

% For this section, use RANSAC to find the best fundamental matrix by
% randomly sample interest points. You would reuse
% estimate_fundamental_matrix() from part 2 of this assignment.

% If you are trying to produce an uncluttered visualization of epipolar
% lines, you may want to return no more than 30 points for either left or
% right images.

function [ Best_Fmatrix, inliers_a, inliers_b] = ransac_fundamental_matrix(matches_a, matches_b)


%%%%%%%%%%%%%%%%
% Your code here
itera = 2000; %% for Mont and Gaudi 500
sigma=0.01; %%for unnorm Dotre: 0.008; for unnorm Mont 0.02
%sigma=0.03;  %%for norm Notre: 0.03; for norm Mont 0.06
sample=8;

L=size(matches_a,1);

count=zeros(itera,1);

sample_a=[];
sample_b=[];
temp_F=[];
Best_Fmatrix=zeros(3,3);
Max_inlier=0;
inlier=[];

for K=1:1:itera    
    % step1 - find sample calcu the F
    cc=randperm(L);
    c=cc(1:sample);
    sample_a=matches_a(c,:);
    sample_b=matches_b(c,:);
    temp_F = estimate_fundamental_matrix(sample_a,sample_b);
    %temp_F = estimate_fundamental_matrix_normal(sample_a,sample_b);
    % step2 - calcu the ri for every point 
    ri=[];
    for m=1:1:L
        xa=[matches_a(m,:)';1];
        xb=[matches_b(m,:)';1];
        %xbpre=temp_F*xa;        
        %dis=(sum((xb(1:2)-xbpre(1:2)).^2))^(1/2);
        dis=sum(abs(xb'*temp_F*xa).^2)^(1/2);
        ri=[ri;dis];
    end
    % step3 - decide inliers
    fin=ri-sigma;
    ind=find(fin<=0);
    Ninlier=length(ind);
    count(K)=Ninlier;
    % step4 - update max and fill in inlier
    if Ninlier>Max_inlier
        Max_inlier=Ninlier;
        Best_Fmatrix=temp_F;
        inlier=ind;
    else
    end    
end
%%%%%%%%%%%%%%%%
inliers_a=matches_a(inlier,:);
inliers_b=matches_b(inlier,:);
% Your ransac loop should contain a call to 'estimate_fundamental_matrix()'
% that you wrote for part II.

%placeholders, you can delete all of this
% Best_Fmatrix = estimate_fundamental_matrix(matches_a(1:10,:), matches_b(1:10,:));
% inliers_a = matches_a(1:30,:);
% inliers_b = matches_b(1:30,:);
end

