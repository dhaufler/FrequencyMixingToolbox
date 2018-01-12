function [obj] = FM_TripletClusterCentroids(obj,threshold)
%Use with output from 'WrappedPhaseTriplets.m'
%   To do: Output should clearly differentiate b/w units of frequency and
%   indices

%init parameters
n_freqs = obj.N_freqs;
freqs = obj.freqs;

%output fields
obj.circR_analysis.XYZ_peaks = cell(8,1);
obj.circR_analysis.IJ_peaks = cell(8,1);
obj.circR_analysis.IJ_data = zeros(n_freqs,n_freqs,8);
obj.circR_analysis.IJ_data_no_control = zeros(n_freqs,n_freqs,8);
obj.circR_analysis.IJ_circR_filtered = [];

%For each of the 8 cluster-types, convert the circR scores onto a 2D grid,
%using 2 axes. The axes are selected to be the 2 components of X, Y, and Z
%with the most dense grid spacing:
IJ_inds = [     2 3; ...1
                2 3; ...2
                2 3; ...3
                2 3; ...4
                1 3; ...5
                1 2; ...6
                2 3; ...7
                2 3; ...8
                1 2; ...9
                1 3; ...10
                2 3; ...11
                1 2; ...12
                2 3]; %13
  
for i=1:13 %for each relationship type
    for j = 1:length(obj.circR_analysis.Test.circR{i})
        I_ind = obj.circR_analysis.Frequency_Indices{i}(j,IJ_inds(i,1));
        J_ind = obj.circR_analysis.Frequency_Indices{i}(j,IJ_inds(i,2));
        
        obj.circR_analysis.IJ_data_no_control(I_ind,J_ind,i)= obj.circR_analysis.Test.circR{i}(j);
        obj.circR_analysis.IJ_data(I_ind,J_ind,i,1) = obj.circR_analysis.Test.circR{i}(j)-max(obj.circR_analysis.Control.circR{i}(:,j));
        obj.circR_analysis.IJ_data(I_ind,J_ind,i,2) = obj.circR_analysis.Test.circRp{i}(j);
    end
end

%%Apply filter to each relationship type
%for i=1:8
%    IJ_circR_filtered(:,:,i) = medfilt2(squeeze(IJ_data(:,:,i,1)),[3 3]);
%end
obj.circR_analysis.IJ_circR_filtered = squeeze(obj.circR_analysis.IJ_data(:,:,:,1));

% get coordinates of maxima of filtered data
for i=1:13
    filt = (fspecial('gaussian', 7,5));
    %p = FastPeakFind(squeeze(IJ_circR_filtered(:,:,i,1)),0.01,filt,1);
    p = FastPeakFind(squeeze(obj.circR_analysis.IJ_circR_filtered(:,:,i,1)),0.01,filt);
    I_ind = p(2:2:end);
    J_ind = p(1:2:end);
    CR_val = zeros(length(I_ind),1);
    for j=1:length(I_ind)
        CR_val(j,1) = obj.circR_analysis.IJ_circR_filtered(round(I_ind(j)),round(J_ind(j)),i,1);
    end
    %keep peaks above threshold
    good_peaks = CR_val>threshold; %mean resultant vector threshold
    obj.circR_analysis.IJ_peaks{i} = [I_ind(good_peaks) J_ind(good_peaks) CR_val(good_peaks)];
end

% %Finally, convert IJ_peaks into X,Y,Z values in frequency units
% 
% k=1;
% if size(obj.circR_analysis.IJ_peaks{k},1)>0
% temp_X = obj.circR_analysis.IJ_peaks{k}(:,1) + obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_Y = obj.circR_analysis.IJ_peaks{k}(:,1);
% temp_Z = obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_CRval = obj.circR_analysis.IJ_peaks{k}(:,3);
% [~,ind_order] = sort(temp_CRval,1,'descend');
% temp = [temp_X temp_Y temp_Z temp_CRval];
% temp2 = temp(ind_order,:);
% XYZ_peaks{k}=interp1(1:n_freqs,freqs,temp2,'linear');
% end
% 
% k=2;
% if size(obj.circR_analysis.IJ_peaks{k},1)>0
% temp_X = obj.circR_analysis.IJ_peaks{k}(:,1) + 2*obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_Y = obj.circR_analysis.IJ_peaks{k}(:,1);
% temp_Z = obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_CRval = obj.circR_analysis.IJ_peaks{k}(:,3);
% [~,ind_order] = sort(temp_CRval,1,'descend');
% temp = [temp_X temp_Y temp_Z temp_CRval];
% temp2 = temp(ind_order,:);
% XYZ_peaks{k}=interp1(1:n_freqs,freqs,temp2,'linear');
% end
% 
% k=3;
% if size(obj.circR_analysis.IJ_peaks{k},1)>0
% temp_X = 2*obj.circR_analysis.IJ_peaks{k}(:,1) + obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_Y = obj.circR_analysis.IJ_peaks{k}(:,1);
% temp_Z = obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_CRval = obj.circR_analysis.IJ_peaks{k}(:,3);
% [~,ind_order] = sort(temp_CRval,1,'descend');
% temp = [temp_X temp_Y temp_Z temp_CRval];
% temp2 = temp(ind_order,:);
% obj.circR_analysis.XYZ_peaks{k}=interp1(1:n_freqs,freqs,temp2,'linear');
% end
% 
% k=4;
% if size(obj.circR_analysis.IJ_peaks{k},1)>0
% temp_X = 2*obj.circR_analysis.IJ_peaks{k}(:,1) - obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_Y = obj.circR_analysis.IJ_peaks{k}(:,1);
% temp_Z = obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_CRval = obj.circR_analysis.IJ_peaks{k}(:,3);
% [~,ind_order] = sort(temp_CRval,1,'descend');
% temp = [temp_X temp_Y temp_Z temp_CRval];
% temp2 = temp(ind_order,:);
% obj.circR_analysis.XYZ_peaks{k}=interp1(1:n_freqs,freqs,temp2,'linear');
% end
% 
% k=5;
% if size(obj.circR_analysis.IJ_peaks{k},1)>0
%     temp_X = obj.circR_analysis.IJ_peaks{k}(:,1);
%     temp_Y = 2*obj.circR_analysis.IJ_peaks{k}(:,1)-2*obj.circR_analysis.IJ_peaks{k}(:,2);
%     temp_Z = obj.circR_analysis.IJ_peaks{k}(:,2);
%     temp_CRval = obj.circR_analysis.IJ_peaks{k}(:,3);
%     [~,ind_order] = sort(temp_CRval,1,'descend');
%    temp = [temp_X temp_Y temp_Z temp_CRval];
% temp2 = temp(ind_order,:);
% obj.circR_analysis.XYZ_peaks{k}=interp1(1:n_freqs,freqs,temp2,'linear');
% end
% 
% k=6;
% if size(obj.circR_analysis.IJ_peaks{k},1)>0
% temp_X = obj.circR_analysis.IJ_peaks{k}(:,1);
% temp_Y = obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_Z = 2*obj.circR_analysis.IJ_peaks{k}(:,1)-2*obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_CRval = obj.circR_analysis.IJ_peaks{k}(:,3);
% [~,ind_order] = sort(temp_CRval,1,'descend');
% temp = [temp_X temp_Y temp_Z temp_CRval];
% temp2 = temp(ind_order,:);
% obj.circR_analysis.XYZ_peaks{k}=interp1(1:n_freqs,freqs,temp2,'linear');
% end
% 
% k=7;
% if size(obj.circR_analysis.IJ_peaks{k},1)>0
% temp_X = 2*obj.circR_analysis.IJ_peaks{k}(:,1)+2*obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_Y = obj.circR_analysis.IJ_peaks{k}(:,1);
% temp_Z = obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_CRval = obj.circR_analysis.IJ_peaks{k}(:,3);
% [~,ind_order] = sort(temp_CRval,1,'descend');
% temp = [temp_X temp_Y temp_Z temp_CRval];
% temp2 = temp(ind_order,:);
% obj.circR_analysis.XYZ_peaks{k}=interp1(1:n_freqs,freqs,temp2,'linear');
% end
% 
% k=8;
% if size(obj.circR_analysis.IJ_peaks{k},1)>0
% temp_X = 2*obj.circR_analysis.IJ_peaks{k}(:,1)-2*obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_Y = obj.circR_analysis.IJ_peaks{k}(:,1);
% temp_Z = obj.circR_analysis.IJ_peaks{k}(:,2);
% temp_CRval = obj.circR_analysis.IJ_peaks{k}(:,3);
% [~,ind_order] = sort(temp_CRval,1,'descend');
% temp = [temp_X temp_Y temp_Z temp_CRval];
% temp2 = temp(ind_order,:);
% obj.circR_analysis.XYZ_peaks{k}=interp1(1:n_freqs,freqs,temp2,'linear');
% end
end

