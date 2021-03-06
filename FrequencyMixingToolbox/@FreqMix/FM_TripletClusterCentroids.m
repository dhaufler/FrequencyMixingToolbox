function [obj] = FM_TripletClusterCentroids(obj,threshold)
%Use with output from 'WrappedPhaseTriplets.m'
%   To do: Output should clearly differentiate b/w units of frequency and
%   indices

%init parameters
n_freqs     = obj.N_freqs;
freqs       = obj.freqs;

%output fields
obj.circR_analysis.XYZ_peaks            = cell(obj.N_distinct_XYZ_relationships,1);
obj.circR_analysis.IJ_peaks             = cell(obj.N_distinct_XYZ_relationships,1);
obj.circR_analysis.IJ_data              = NaN(n_freqs,n_freqs,obj.N_distinct_XYZ_relationships);
obj.circR_analysis.IJ_Surrogate         = NaN(n_freqs,n_freqs,obj.N_distinct_XYZ_relationships);
obj.circR_analysis.IJ_data_no_control   = NaN(n_freqs,n_freqs,obj.N_distinct_XYZ_relationships);
obj.circR_analysis.IJ_circR_filtered    = [];

%For each of the 8 cluster-types, convert the circR scores onto a 2D grid,
%using 2 axes. The axes are selected to be the 2 components of X, Y, and Z
%with the most dense grid spacing:
% IJ_inds = [     2 3; ...1
%                 2 3; ...2
%                 2 3; ...3
%                 2 3; ...4
%                 1 3; ...5
%                 1 2; ...6
%                 2 3; ...7
%                 2 3; ...8
%                 1 2; ...9
%                 1 3; ...10
%                 2 3; ...11
%                 1 2; ...12
%                 2 3]; %13
  
for i=1:obj.N_distinct_XYZ_relationships %for each relationship type
    for j = 1:length(obj.circR_analysis.Test.circR{i})
        I_ind = obj.circR_analysis.Frequency_Indices{i}(j,2);
        J_ind = obj.circR_analysis.Frequency_Indices{i}(j,3);
        
        obj.circR_analysis.IJ_data_no_control(I_ind,J_ind,i)= obj.circR_analysis.Test.circR{i}(j);
        obj.circR_analysis.IJ_data(I_ind,J_ind,i,1) = obj.circR_analysis.Test.circR{i}(j)-obj.circR_analysis.Surrogate.circR{i}(j);
        obj.circR_analysis.IJ_data(I_ind,J_ind,i,2) = obj.circR_analysis.Test.circRp{i}(j);
    end
end

% %Apply filter to each relationship type
% for i=1:obj.N_distinct_XYZ_relationships
%    IJ_circR_filtered(:,:,i) = medfilt2(squeeze(IJ_data(:,:,i,1)),[3 3]);
% end
obj.circR_analysis.IJ_circR_filtered = squeeze(obj.circR_analysis.IJ_data(:,:,:,1));

% get coordinates of maxima of filtered data
% for i=1:obj.N_distinct_XYZ_relationships
%     filt = (fspecial('gaussian', 7,5));
%     %p = FastPeakFind(squeeze(IJ_circR_filtered(:,:,i,1)),0.01,filt,1);
%     p = FastPeakFind(squeeze(obj.circR_analysis.IJ_circR_filtered(:,:,i,1)),0.01,filt);
%     I_ind = p(2:2:end);
%     J_ind = p(1:2:end);
%     CR_val = zeros(length(I_ind),1);
%     for j=1:length(I_ind)
%         CR_val(j,1) = obj.circR_analysis.IJ_circR_filtered(round(I_ind(j)),round(J_ind(j)),i,1);
%     end
%     %keep peaks above threshold
%     good_peaks = CR_val>threshold; %mean resultant vector threshold
%     obj.circR_analysis.IJ_peaks{i} = [I_ind(good_peaks) J_ind(good_peaks) CR_val(good_peaks)];
% end

end

