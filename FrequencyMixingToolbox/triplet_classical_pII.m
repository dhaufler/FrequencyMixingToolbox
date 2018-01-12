function [H] = triplet_classical_pII(X_)
%UNTITLED Summary of this function goes here
%   X is (N x 3) complex values, where N is the number of (independent) samples. 

%parameters and init
X = angle(X_);
N_voxel_edges = 5;
N_phase_samples = size(X,1);
H_data = [];
BUB_coeffs = []; %for BUB entropy estimator
B = []; % for histndim


%get the BUB coefficients
[BUB_coeffs.a1,~]=modBUBfunc(N_phase_samples,N_voxel_edges,11,0);
[BUB_coeffs.a2,~]=modBUBfunc(N_phase_samples,N_voxel_edges.^2,11,0);
[BUB_coeffs.a3,~]=modBUBfunc(N_phase_samples,N_voxel_edges.^3,11,0);

%create the B parameter matrices for histndim
%for each of the 1,2, and 3 dimensional entropies, we require a parameter matrix B for the function histndim.m
B.B_1 = zeros(3,1);
B.B_1(1,:) = ones(1,1)*N_voxel_edges;
B.B_1(2,:) = ones(1,1)*(-pi);
B.B_1(3,:) = ones(1,1)*(pi);

B.B_12 = zeros(3,2);
B.B_12(1,:) = ones(1,2)*N_voxel_edges;
B.B_12(2,:) = ones(1,2)*(-pi);
B.B_12(3,:) = ones(1,2)*(pi);

B.B_123 = zeros(3,3);
B.B_123(1,:) = ones(1,3)*N_voxel_edges;
B.B_123(2,:) = ones(1,3)*(-pi);
B.B_123(3,:) = ones(1,3)*(pi);

% Compute pII (this code is copied from 'myInteractionInformation'
%get binned histogram counts
[N_X,~]=histndim(X(:,1),B.B_1);          
[N_Y,~]=histndim(X(:,2),B.B_1);
[N_Z,~]=histndim(X(:,3),B.B_1);
[N_XY,~]=histndim(X(:,[1 2]),B.B_12);
[N_YZ,~]=histndim(X(:,[2 3]),B.B_12);
[N_XZ,~]=histndim(X(:,[1 3]),B.B_12);
[N_XYZ,~]=histndim(X(:,[1 2 3]),B.B_123);


%compute the BUB entropy estimate
H.X = sum(BUB_coeffs.a1(N_X(:)+1));
H.Y = sum(BUB_coeffs.a1(N_Y(:)+1));
H.Z = sum(BUB_coeffs.a1(N_Z(:)+1));
H.XY = sum(BUB_coeffs.a2(N_XY(:)+1));
H.YZ = sum(BUB_coeffs.a2(N_YZ(:)+1));
H.XZ = sum(BUB_coeffs.a2(N_XZ(:)+1));
H.XYZ = sum(BUB_coeffs.a3(N_XYZ(:)+1));

%compute interaction information
H.II = H.XY + H.YZ + H.XZ - H.X - H.Y - H.Z - H.XYZ;

end

