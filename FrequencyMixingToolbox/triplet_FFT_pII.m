function [H] = triplet_FFT_pII(Triplet_Spectral_Data)
%UNTITLED Summary of this function goes here
% Input:    'Triplet_Spectral_Data' is (N x 3) complex values, where N is the 
%           number of (independent) samples. 
% Output:   'H' is a struct with fields 'pII', a scalar, and 'S_rel', an
%           8x4 vector of the 8 frequency/phase relationships (first three columns) together with
%           the corresponding FFT value, and 'S_pw',
%           a 3x1 vector of the largest FFT values in each of the three
%           "pairwise-only dimensions 

%parameters and init
X = angle(Triplet_Spectral_Data);
N_voxel_edges = 5;
H = [];
H.pII = NaN;
H.S_rel = NaN(8,4);
H.S_pw = NaN(3,1);

%get binned histogram counts
N_XYZ = histcn(X, linspace(-pi,pi,N_voxel_edges+1), linspace(-pi,pi,N_voxel_edges+1), linspace(-pi,pi,N_voxel_edges+1));

% Compute 3D FFT and center the result. Frequency indices for S are as
% follows:
% In X dimension: 0 1 2
% In Y dimension: -2 -1 0 1 2
% In Z dimension: -2 -1 0 1 2
S__ = fftn((N_XYZ-mean(N_XYZ))/sum(N_XYZ(:)),size(N_XYZ));
S_ = fftshift(S__);
S = abs(S_(3:5,:,:)).^2;

%get each of the "pairwise" planes, where the x, y, or z frequency is 0;
S_x0 = squeeze(S(1,:,:));
S_y0 = squeeze(S(:,3,:));
S_z0 = squeeze(S(:,:,3));

%get the max in each of the above planes
H.S_pw(1) = max(S_x0(:));
H.S_pw(2) = max(S_y0(:));
H.S_pw(3) = max(S_z0(:));

%get output
H.pII = max(S(:)) - sum(H.S_pw);

%get each relationship type. The order of the 8 types is based on
%'Triplet_Table', saved in 'Triplet_Table_standard.mat'
H.S_rel(1,:) = [1, -1, -1, S(2,2,2)];
H.S_rel(2,:) = [1, -1, -2, S(2,2,1)];
H.S_rel(3,:) = [1, -2,  1, S(2,1,4)];
H.S_rel(4,:) = [1, -2, -1, S(2,1,2)];
H.S_rel(5,:) = [2, -1, -2, S(3,2,1)];
H.S_rel(6,:) = [2, -2, -1, S(3,1,2)];
H.S_rel(7,:) = [1, -2, -2, S(2,1,1)];
H.S_rel(8,:) = [1, -2,  2, S(2,1,5)]; 

end