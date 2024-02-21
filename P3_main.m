%> ------------------------------------------------------------------------
%> ENGN2560: Computer Vision
%>    Lab04: Absolute Camera Pose and Visual Odometry
%> Problem3: Visual Odometry Part II
%> ------------------------------------------------------------------------
clc; clear all; close all;
rng(0);

%> Load camera intrinsic matrix (borrowed from Problem2)
load('MyData/Problem2/IntrinsicMatrix.mat');

%> Load ground truth poses (borrowed from Problem2)
load('data/Problem2/GT_Poses.mat');	%> GT_Poses

%> Read all images in the sequence.
%> Use imread(Image_Sequence(i).name); to read image i
mfiledir = fileparts(mfilename('fullpath'));
Image_Sequence = dir([mfiledir, '/data/Problem2/fr2_desk/*.png']);

%> Parameters passed to RANSAC
PARAMS.INLIER_THRESH                          = 2;      %> 2 pixels
PARAMS.RANSAC_ITERATIONS                      = 1500;   %> Total number of RANSAC iterations
PARAMS.NUM_OF_FRAMES_FROM_LAST_KF             = 20;
PARAMS.RATIO_OF_COVISIBLE_POINTS_FROM_LAST_KF = 0.5;

%> =================================================================
%> TODO: Implement a Visual Odometry Pipeline Reducing Motion Drift
%> =================================================================
nImgs = length(Image_Sequence);

Img1 = imread(Image_Sequence(1).name);
Img2 = imread(Image_Sequence(2).name);

% code from part 1



% computing scale

Rgi = GT_Poses(:,1:3,i);
Tgi = GT_Poses(:,4,i);
cgi = -Rgi' * Tgi;

Rgim1 = GT_Poses(:,1:3,i-1);
Tgim1 = GT_Poses(:,4,i-1);
cgim1 = -Rgim1' * Tgim1;

s = norm(cgi - cgim1);

% rest of the images

prev_kf_triangulated;
prev_kf_points3D;
prev_kf_idx = 2;

kf_idxs = zeros(1, nImgs);
kf_idxs(1:2) = [1,1];

estimatedPoses = zeros(3, 4, nImgs);
estimatedPoses(:,:,1) = [eye(3), zeros(3,1)];

for i=3:nImgs
    currImg = imread(Image_Sequence(i).name);
    
    % compute absolute camera pose
    [AbsR, AbsT] = Ransac4AbsPose(PARAMS, points3D, points2D, K);
    
    Np = 0;
    Nq = 0;
    Nf = Np + Nq;
    covisRatio = Nq / Nf;
    
    if ((i - prev_kf_idx) > PARAMS.NUM_OF_FRAMES_FROM_LAST_KF) || ...
        (covisRatio < PARAMS.RATIO_OF_COVISIBLE_POINTS_FROM_LAST_KF)
        
        
        prev_kf = currImg;
        prev_kf_idx = i;
        
        kf_idxs(i) = 1;
    end
    
    estimatedPoses(:,:,i) = [AbsR, s*AbsT];
end

Visualize_Trajectory(GT_Poses, estimatedPoses, find(kf_idxs == 1));

% compute RMSE for R and T
r_sum = 0;
t_sum = 0;
for i=1:nImgs
    R_gt = GT_Poses(:, 1:3, i);
    T_gt = GT_Poses(:, 4, i);
    
    R = estimatedPoses(:, 1:3, i);
    T = estimatedPoses(:, 4, i);
    
    [deltaR, deltaT] = getRT_Error(R, T, R_gt, T_gt);
    r_sum = r_sum + (deltaR^2);
    t_sum = t_sum + (deltaT^2);
end

rmse_R = sqrt(r_sum / nImgs);
rmse_T = sqrt(t_sum / nImgs);

disp(rmse_R);
disp(rmse_T);
