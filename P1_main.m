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

[R, T, E, matched1, matched2, matches, inlierIdxs] = estimateRelativePose(PARAMS, Img1, Img2, K);

R1 = eye(3);
T1 = zeros(3,1);
R2 = R;
T2 = T;

nPoints = size(inlierIdxs, 2);

pointsSymmedian = zeros(3, nPoints);

for i=1:nPoints
    gamma1 = [matched1(:, i); 1];
    gamma2 = [matched2(:, i); 1];
    
    worldPoint = symmedianTriangulation(gamma1, gamma2, R1, T1, R2, T2);
    pointsSymmedian(:, i) = worldPoint;
end

img2_bw = single(rgb2gray(Img2));
img3_bw = single(rgb2gray(Img3));

[f2, d2] = vl_sift(img2_bw);
[f3, d3] = vl_sift(img3_bw);

[matches23, scores] = vl_ubcmatch(d2, d3); % TODO: DEBUG THIS PART

matches12_f2_idxs = matches(2, inlierIdxs); % features that are inliers

inMatches23 = intersect(matches23(1,:), matches12_f2_idxs); % inlier features and present in 3
inMatched2 = ismember(matches12_f2_idxs, inMatches23); % whether inlier features are present in 3

points2D = f3(1:2, inMatches23);
points3D = pointsSymmedian(:, inMatched2);

[AbsR, AbsT] = Ransac4AbsPose(PARAMS, points3D, points2D, K);

[deltaR, deltaT] = getRT_Error(AbsR, AbsT, R_gt, T_gt);

disp(deltaR);
disp(deltaT);

% FUNCTIONS ---------------------------------------------------------------

function [deltaR, deltaT] = getRT_Error(R, T, Rg, Tg)
    deltaR = acos((trace(Rg' * R) - 1) / 2);

    Tg = Tg / norm(Tg);
    T = T / norm(T);

    deltaT = abs(dot(Tg, T) - 1);
end

function S = getSkewSymmetric(v)
    S = [0, -v(3), v(2);
     v(3), 0, -v(1);
     -v(2), v(1), 0];
end

function [AbsR, AbsT] = Ransac4AbsPose(PARAMS, Points3D, Points2D, K)
    maxInliers = 0;
    
    Points3DHomogenous = [Points3D; ones(1, size(Points3D, 2))]; % TODO: CAN I DO THIS?
    
    for i=1:PARAMS.RANSAC_ITERATIONS
        sample_col_idxs = randperm(size(Points3D, 2), 3);
        sample2D = [Points2D(:, sample_col_idxs); ones(1, 3)];
        sample3D = Points3D(:, sample_col_idxs);
        
%         disp(samples);
        
        [Rs, Ts] = P3P_LambdaTwist(sample2D, sample3D);
        
        num_candidates = size(Ts, 2);
        for j=1:num_candidates
            R_test = Rs(:, :, j);
            T_test = Ts(:, j);
            RT_test = [R_test, T_test];
       
            test_M = K * RT_test;
            reprojected = test_M * Points3DHomogenous;
%             reprojected(1,:) = reprojected(1,:) ./ reprojected(3,:);
%             reprojected(2,:) = reprojected(2,:) ./ reprojected(3,:);
            
%             disp(reprojected);

            distances = sqrt(sum((reprojected(1:2,:) - Points2D).^2));
            
%             disp(distances);
            
            inliers = distances < PARAMS.INLIER_THRESH;
            
            inlierCount = sum(inliers);
            
            if inlierCount > maxInliers
%                 disp(essential_test);
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
% 
% function distance = pointToPlaneDistance(point, PARAMS)
%     normOfPlaneNorm = norm(PARAMS.GT_PLANE_NORMAL_VECTOR);
%     distance = abs(dot(PARAMS.GT_PLANE_NORMAL_VECTOR, point) + PARAMS.GT_PLANE_DISPLACEMENT) / normOfPlaneNorm;
% end

function [R, T, E] = deriveTransformations(E_final, matched1, matched2)
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

    disp(det(R1));
    disp(det(R2));
    
    candidateRTs = {R1, T_plus, R2, T_plus, R1, T_minus, R2, T_minus};

    mostPos = 0;
    bestR = [];
    bestT = [];
    for i=1:2:length(candidateRTs)
        curr_R = candidateRTs{i};
        curr_T = candidateRTs{i+1};
        posCount = 0;
        for j=1:size(matched1, 2)
            feature_1_homogenous = [matched1(:, j); 1];
            feature_2_homogenous = [matched2(:, j); 1];

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

function [R, T, E, matched1, matched2, matches, inlierIdxs] = estimateRelativePose(PARAMS, Img1, Img2, K)

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
    [~, sortedIndices] = sort(scores, 'descend');
    sortedMatches = matches(:, sortedIndices);
    n = floor(0.8 * size(sortedMatches, 2));
    top_n = sortedMatches(:, 1:n);
    
%     top_n = matches(:, 1:n);

    matched1 = f1(1:2, top_n(1, :));
    matched2 = f2(1:2, top_n(2, :));

    [E_final, inlierIdxs] = Ransac4Essential(PARAMS, matched1, matched2, K);
    
    matched1 = matched1(:, inlierIdxs);
    matched2 = matched2(:, inlierIdxs);

    [R,T,E] = deriveTransformations(E_final, matched1, matched2);

end

% RANSAC ------------------------------------------------------------------

function [E, inlier_Idx] = Ransac4Essential(PARAMS, gamma1, gamma2, K)
    maxInliers = 0;
    
    gamma1_homogeneous = [gamma1; ones(1, size(gamma1, 2))];
    gamma2_homogeneous = [gamma2; ones(1, size(gamma2, 2))];
    
    for i=1:PARAMS.RANSAC_ITERATIONS
        sample_col_idxs = randperm(size(gamma1, 2), 5);
        samples = ones(5, 3, 2);
        samples(:, 1:2, 1) = gamma1(:, sample_col_idxs)';
        samples(:, 1:2, 2) = gamma2(:, sample_col_idxs)';
        
%         disp(samples);
        
        essential_candidates = fivePointAlgorithmSelf(samples);
        
        for j=1:length(essential_candidates)
            essential_test = essential_candidates{j};
            F = transpose(inv(K)) * essential_test / K;
            
            epipolar_lines = F * gamma1_homogeneous;
            distances = abs(sum(epipolar_lines .* gamma2_homogeneous, 1)) ./ sqrt(sum(epipolar_lines(1:2, :).^2, 1));
            
            inliers = distances < PARAMS.INLIER_THRESH;
            
            inlierCount = sum(inliers);
            
            if inlierCount > maxInliers
%                 disp(essential_test);
                maxInliers = inlierCount;                
                inlier_Idx = find(inliers == 1);
                E = essential_test;
            end
        end
    end
end

