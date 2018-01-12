function [obj] = FM_FFT_pII(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Init
phase_bins = linspace(-pi,pi,6);
obj.fft_pII_analysis.Test = [];
obj.fft_pII_analysis.Surrogate = [];

% FFT indices (for shifted matrix) for 8 relationship types as ordered in
% Triplet_Table
FFT_inds = [4 2 2;...
            4 2 1;...
            4 1 4;...
            4 1 2;...
            5 2 1;...
            5 1 2;...
            4 1 1;...
            4 1 5];
        
%for triplets lying on the planes of interest, compute 3D fft.

%for shuffling
X_shuffled_inds = randsample(size(obj.data,1),size(obj.data,1));
for xxx = 1:obj.N_distinct_XYZ_relationships %for each F-P relationship
    %Data:
    X_data = angle(obj.data(:,obj.Frequency_Indices{xxx}(:,1)));
    Y_data = angle(obj.data(:,obj.Frequency_Indices{xxx}(:,2)));
    Z_data = angle(obj.data(:,obj.Frequency_Indices{xxx}(:,3)));
    
    %shuffle 
    %X_data = X_data(X_shuffled_inds,:);
    
    %Surrogate Data:
   % X_data_S = angle(obj.S_data(:,obj.Frequency_Indices{xxx}(:,1)));
    %Y_data_S = angle(obj.S_data(:,obj.Frequency_Indices{xxx}(:,2)));
   % Z_data_S = angle(obj.S_data(:,obj.Frequency_Indices{xxx}(:,3)));
    
    %init output structures
    obj.fft_pII_analysis.Test.S{xxx,1} = zeros(size(obj.Frequency_Indices{xxx},1),5,5,5);
    obj.fft_pII_analysis.Surrogate.S{xxx,1} = zeros(size(obj.Frequency_Indices{xxx},1),5,5,5);
    temp_out_S          = zeros(size(obj.Frequency_Indices{xxx},1),5,5,5);
    temp_out_S_surr     = zeros(size(obj.Frequency_Indices{xxx},1),5,5,5);
    temp_out_pII        = zeros(size(obj.Frequency_Indices{xxx},1),2);
    temp_out_pII_surr   = zeros(size(obj.Frequency_Indices{xxx},1),2);
    temp_out_pIIs        = zeros(size(obj.Frequency_Indices{xxx},1),2); %pII-specific--compute pII with the expected spectral index
    temp_out_pIIs_surr   = zeros(size(obj.Frequency_Indices{xxx},1),2);
    for ij = 1:size(obj.Frequency_Indices{xxx},1) %for each index
        %for data
        temp_binned_data1 = histcn([X_data(:,ij),Y_data(:,ij),Z_data(:,ij)],...
            phase_bins, phase_bins, phase_bins);
        temp_binned_data2 = (temp_binned_data1 - mean(temp_binned_data1(:)))/sum(temp_binned_data1(:));
        
        S__ = fftn(temp_binned_data2,[5 5 5]);
        S_ = fftshift(S__);
        S = abs(S_).^2;
        S_xy = squeeze(S(:,:,3));
        S_xz = squeeze(S(:,3,:));
        S_yz = squeeze(S(3,:,:));
        S_specific = squeeze(S(FFT_inds(xxx,1),FFT_inds(xxx,2),FFT_inds(xxx,3)));       
        temp_out_S(ij,:,:,:) = S;
        [max_S,temp_out_pII(ij,2)] = max(S(:));
        temp_out_pII(ij,1) = max_S - min(sqrt([max(S_xy(:)), max(S_xz(:)), max(S_yz(:))]));
        temp_out_pIIs(ij,1) = S_specific - min(sqrt([max(S_xy(:)), max(S_xz(:)), max(S_yz(:))]));      
        
%         %for surrogate data
%         temp_binned_data1_S = histcn([X_data_S(:,ij),Y_data_S(:,ij),Z_data_S(:,ij)],...
%             phase_bins, phase_bins, phase_bins);
%         temp_binned_data2_S = (temp_binned_data1_S - mean(temp_binned_data1_S(:)))/sum(temp_binned_data1_S(:));
%         
%         S__ = fftn(temp_binned_data2_S,[5 5 5]);
%         S_ = fftshift(S__);
%         S = abs(S_).^2;
%         S_xy = squeeze(S(:,:,1));
%         S_xz = squeeze(S(:,1,:));
%         S_yz = squeeze(S(1,:,:));
%         S_specific = squeeze(S(FFT_inds(xxx,1),FFT_inds(xxx,2),FFT_inds(xxx,3)));
%         temp_out_S_surr(ij,:,:,:) = S;
%         [max_S,temp_out_pII_surr(ij,2)] = max(S(:));
%         temp_out_pII_surr(ij,1) = max_S - min(sqrt([max(S_xy(:)), max(S_xz(:)), max(S_yz(:))]));
%         temp_out_pIIs_surr(ij,1) = S_specific - min(sqrt([max(S_xy(:)), max(S_xz(:)), max(S_yz(:))]));    
        
    end
    obj.fft_pII_analysis.Test.S{xxx,1} = temp_out_S;
    obj.fft_pII_analysis.Test.pII{xxx,1} = temp_out_pII;
    obj.fft_pII_analysis.Test.pIIs{xxx,1} = temp_out_pIIs;
%     obj.fft_pII_analysis.Surrogate.S{xxx,1} = temp_out_S_surr;
%     obj.fft_pII_analysis.Surrogate.pII{xxx,1} = temp_out_pII_surr;   
%     obj.fft_pII_analysis.Surrogate.pIIs{xxx,1} = temp_out_pIIs_surr;  
end     
end

