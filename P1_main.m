%> ------------------------------------------------------------------------
%> ENGN2560: Computer Vision
%>    Lab04: Absolute Camera Pose and Visual Odometry
%> Problem1: Absolute Pose Estimation
%> ------------------------------------------------------------------------
clc; clear all; close all;
rng(0);

%> Read triplet image pair
Img1 = imread('data/Problem1/Img1.png');
Img2 = imread('data/Problem1/Img2.png');
Img3 = imread('data/Problem1/Img3.png');

%> Load camera intrinsic matrix
load('data/Problem1/IntrinsicMatrix.mat');

%> Load ground truth pose for the third image
load('data/Problem1/GT_R.mat');		%> R_gt
load('data/Problem1/GT_T.mat');		%> T_gt

%> Parameters passed to RANSAC
PARAMS.INLIER_THRESH                 = 2;      %> 2 pixels
PARAMS.RANSAC_ITERATIONS             = 1000;   %> Total number of RANSAC iterations

%> =========================================================
%> TODO: Estimate an Absolute Pose under a RANSAC scheme
%> =========================================================

[R, T, E, matched1, matched2, d2, inlierMatches12, inlierIdxs] = estimateRelativePose(PARAMS, Img1, Img2, K);

R1 = eye(3);
T1 = zeros(3,1);
R2 = R;
T2 = T;

nPoints = size(inlierIdxs, 2);

pointsSymmedian = zeros(3, nPoints);

for i=1:nPoints
    gamma1 = K \ [matched1(:, i); 1];
    gamma2 = K \ [matched2(:, i); 1];
    
    worldPoint = symmedianTriangulation(gamma1, gamma2, R1, T1, R2, T2);
    pointsSymmedian(:, i) = worldPoint;
end

img3_bw = single(rgb2gray(Img3));

[f3, d3] = vl_sift(img3_bw);

[matches23, scores] = vl_ubcmatch(d2, d3); % TODO: DEBUG THIS PART

[matches12_unique, ~, uniqueIdxs] = unique(inlierMatches12(2,:), 'stable'); % features that are inliers

idxs_of_unique = true(1, size(uniqueIdxs,1));
next = 1;
for i=1:size(uniqueIdxs,1)
    if uniqueIdxs(i) ~= next
        idxs_of_unique(i) = false;
    else
        next = next + 1;
    end
end

matches12_unique = inlierMatches12(2, idxs_of_unique);

matches_from_12_in_23_that_are_inliers = intersect(matches12_unique, matches23(1,:));
matches_from_23_that_are_in_12 = ismember(matches23(1,:), matches_from_12_in_23_that_are_inliers);

inMatched2 = matches23(1, matches_from_23_that_are_in_12);
inMatched3 = matches23(2, matches_from_23_that_are_in_12);

indices_from_12_that_are_inliers_in_23 = ismember(matches12_unique, matches_from_12_in_23_that_are_inliers);
matches_from_12_that_are_inliers_in_23 = matches12_unique(indices_from_12_that_are_inliers_in_23);

[sortedMatches, sortedIndices] = sort(matches_from_12_that_are_inliers_in_23, 'ascend');
% [~, sortedMatched2Indices] = sort(inMatched2, 'ascend');

% get features that are in 12 and 23

% inMatched3 = inMatched3(sortedMatched2Indices);

points2D = f3(1:2, inMatched3);

pointsSymmedian = pointsSymmedian(:, idxs_of_unique);
points3D = pointsSymmedian(:, indices_from_12_that_are_inliers_in_23);
points3D = points3D(:, sortedIndices);

% Sanity Check ------------------------------------------------------------

% CHECKING MATCHED 2 ------------------------------------------------------

% matched2 = matched2(:, idxs_of_unique);
% matched2 = matched2(:, indices_from_12_that_are_inliers_in_23);
% matched2 = matched2(:, sortedIndices);
% 
% Points3DHomogenous = [points3D; ones(1, size(points3D, 2))];
% % points2DMetric = K \ [points2D; ones(1, size(points2D, 2))];
% 
% RT_test = [R2, T2];
% 
% test_M = K * RT_test;
% reprojected = test_M * Points3DHomogenous;
% reprojected(1,:) = reprojected(1,:) ./ reprojected(3,:);
% reprojected(2,:) = reprojected(2,:) ./ reprojected(3,:);
% 
% distances = sqrt(sum((reprojected(1:2,:) - matched2).^2));
% 
% disp(mean(distances));

% CHECKING MATCHED 3 ------------------------------------------------------

% RT_test = [R_gt, T_gt];
% 
% test_M = K * RT_test;
% reprojected = test_M * Points3DHomogenous;
% reprojected(1,:) = reprojected(1,:) ./ reprojected(3,:);
% reprojected(2,:) = reprojected(2,:) ./ reprojected(3,:);
% 
% distances = sqrt(sum((reprojected(1:2,:) - points2D).^2));
% 
% disp(mean(distances));

% -------------------------------------------------------------------------

[AbsR, AbsT] = Ransac4AbsPose(PARAMS, points3D, points2D, K);

[deltaR, deltaT] = getRT_Error(AbsR, AbsT, R_gt, T_gt);

disp(deltaR);
disp(deltaT);

% Propagating accuracy ----------------------------------------------------

[R3, T3, ~, ~, ~, ~, ~, ~] = estimateRelativePose(PARAMS, Img2, Img3, K);
AbsR3 = R3 * R2;
AbsT3 = R3 * T2 + T3;

[deltaR, deltaT] = getRT_Error(AbsR3, AbsT3, R_gt, T_gt);

disp(deltaR);
disp(deltaT);

% FUNCTIONS ---------------------------------------------------------------

function [deltaR, deltaT] = getRT_Error(R, T, Rg, Tg)
    deltaR = acos((trace(Rg' * R) - 1) / 2);

    Tg = Tg / norm(Tg);
    T = T / norm(T);

    deltaT = abs(dot(Tg, T) - 1);
end

function [AbsR, AbsT] = Ransac4AbsPose(PARAMS, Points3D, Points2D, K)
    maxInliers = 0;

    Points2DHomogenous = [Points2D; ones(1, size(Points2D, 2))];

    Points2D_meters = K \ Points2DHomogenous;
    
    Points3DHomogenous = [Points3D; ones(1, size(Points3D, 2))]; 
    
    for i=1:PARAMS.RANSAC_ITERATIONS
        sample_col_idxs = randperm(size(Points3D, 2), 3);
        sample2D = Points2D_meters(:, sample_col_idxs);
        sample3D = Points3D(:, sample_col_idxs);
                
        [Rs, Ts] = P3P_LambdaTwist(sample2D, sample3D);
        
        num_candidates = size(Ts, 2);
        for j=1:num_candidates
            R_test = Rs(:, :, j);
            T_test = Ts(:, j);
            RT_test = [R_test, T_test];
            % RT_test = [R_gt, T_gt];
       
            test_M = K * RT_test;
            reprojected = test_M * Points3DHomogenous;
            reprojected(1,:) = reprojected(1,:) ./ reprojected(3,:);
            reprojected(2,:) = reprojected(2,:) ./ reprojected(3,:);
            
            distances = sqrt(sum((reprojected(1:2,:) - Points2D).^2));
            % dist = sum(sqrt(sum((reprojected(1:2,:) - points2D').^2, 2)));
            inliers = distances < PARAMS.INLIER_THRESH;
            
            inlierCount = sum(inliers);
            
            if inlierCount > maxInliers
                maxInliers = inlierCount;  
                AbsR = R_test;
                AbsT = T_test;
            end
        end
    end
end

function worldPoint = symmedianTriangulation(gamma1, gamma2, R1, T1, R2, T2)
    c1 = -R1' * T1;
    gamma1NormSq = norm(gamma1) ^ 2;
    inner1 = eye(3) - ((R1' * gamma1 * (gamma1' * R1)) / gamma1NormSq);
    
    c2 = -R2' * T2;
    gamma2NormSq = norm(gamma2) ^ 2;
    inner2 = eye(3) - ((R2' * gamma2 * (gamma2' * R2)) / gamma2NormSq);
    
    worldPoint = (inner1 + inner2) \ ((inner1 * c1) + (inner2 * c2));    
end

function S = getSkewSymmetric(v)
    S = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end

function [R, T, E] = deriveTransformations(E_final, matched1, matched2, K)
    [U, ~, V] = svd(E_final);
    W = [0, -1, 0; 1, 0, 0; 0, 0, 1];
    R1 = U * W * transpose(V);
    R2 = U * transpose(W) * transpose(V);
    
    if (det(R1) < 0 || det(R2) < 0)
        [U, ~, V] = svd(-E_final);
        W = [0, -1, 0; 1, 0, 0; 0, 0, 1];
        R1 = U * W * transpose(V);
        R2 = U * transpose(W) * transpose(V);
    end

    T_plus = U(:, 3);
    T_minus = -U(:, 3);

    % disp(det(R1));
    % disp(det(R2));
    
    candidateRTs = {R1, T_plus, R2, T_plus, R1, T_minus, R2, T_minus};

    mostPos = 0;
    bestR = [];
    bestT = [];
    
    for i=1:2:length(candidateRTs)
        curr_R = candidateRTs{i};
        curr_T = candidateRTs{i+1};

        % disp(curr_R);
        % disp(curr_T);

        posCount = 0;
        for j=1:size(matched1, 2)
            feature_1_homogenous = K \ [matched1(:, j); 1];
            feature_2_homogenous = K \ [matched2(:, j); 1];

            Rg1 = curr_R * feature_1_homogenous;
            A = [-Rg1, feature_2_homogenous];

            res = pinv(A) * curr_T;
            rho1 = res(1);
            rho2 = res(2);

            if (rho1 > 0) && (rho2 > 0)
                posCount = posCount + 1;
            end
        end
        if posCount > mostPos
           mostPos = posCount;
           bestR = curr_R;
           bestT = curr_T;
        end
    end

    bestTx = getSkewSymmetric(bestT); % [0, -bestT(3), bestT(2); bestT(3), 0, -bestT(1); -bestT(2), bestT(1), 0];
    test_E = bestTx * bestR;

    R = bestR;
    T = bestT;
    E = test_E;
end

function [R, T, E, matched1, matched2, d2, top_n, inlierIdxs] = estimateRelativePose(PARAMS, Img1, Img2, K)
    
    % converting to black and white
    img1_bw = single(rgb2gray(Img1));
    img2_bw = single(rgb2gray(Img2));

    % using SFIT
    [f1,d1] = vl_sift(img1_bw);
    [f2,d2] = vl_sift(img2_bw);

    % Matching
    [matches, scores] = vl_ubcmatch(d1, d2);
    % disp(matches);

    % creating rank order list 
    [sortedScores, sortedIndices] = sort(scores, 'ascend');
    sortedMatches = matches(:, sortedIndices);
    n = floor(0.8 * size(sortedMatches, 2));
    top_n = sortedMatches(:, 1:n);
    
    matched1 = f1(1:2, top_n(1, :));
    matched2 = f2(1:2, top_n(2, :));

    [E_final, inlierIdxs] = Ransac4Essential(PARAMS, matched1, matched2, K);
    
    top_n = top_n(:, inlierIdxs);
    matched1 = matched1(:, inlierIdxs);
    matched2 = matched2(:, inlierIdxs);

    [R,T,E] = deriveTransformations(E_final, matched1, matched2, K);

end

% RANSAC ------------------------------------------------------------------

function [E, inlier_Idx] = Ransac4Essential(PARAMS, gamma1, gamma2, K)
    maxInliers = 0;
    
    gamma1_homogeneous = [gamma1; ones(1, size(gamma1, 2))];
    gamma2_homogeneous = [gamma2; ones(1, size(gamma2, 2))];

    gamma1ToMeters = K \ gamma1_homogeneous;
    gamma2ToMeters = K \ gamma2_homogeneous;
    
    for i=1:PARAMS.RANSAC_ITERATIONS
        sample_col_idxs = randperm(size(gamma1, 2), 5);
        samples = ones(5, 3, 2);
        samples(:, 1:2, 1) = gamma1ToMeters(1:2, sample_col_idxs)';
        samples(:, 1:2, 2) = gamma2ToMeters(1:2, sample_col_idxs)';
                
        essential_candidates = fivePointAlgorithmSelf(samples);
        
        for j=1:length(essential_candidates)

            essential_test = essential_candidates{j};
            F = transpose(inv(K)) * essential_test / K;

            epipolar_lines = F * gamma1_homogeneous;
            distances = abs(sum(epipolar_lines .* gamma2_homogeneous, 1)) ./ sqrt(sum(epipolar_lines(1:2, :).^2, 1));

            inliers = distances < PARAMS.INLIER_THRESH;
            
            inlierCount = sum(inliers);
            
            if inlierCount > maxInliers
                maxInliers = inlierCount;                
                inlier_Idx = find(inliers == 1);
                E = essential_test;
            end
        end
    end
end

