function [obj] = FM_R1_R2(obj)
% modified. Original: [R1_R2_freqs,dataOut] = R1_R2_by_triplet_ID(XYZ_freqs,data,triplet_ID,R_max)
%use with 'circR_instead_of_pII'
% Converts a set of frequency indices from XYZ triplets, where X>Y>Z, to
% R1 vs. R2 coordinates, based on one of the 20 triplet IDs.
%   Inputs: 
%   XYZ_freqs should be N-by-3, for N observations
%   data should be N-by-1
%   triplet_ID should be an integer from 1:20
%   R_max: Only consider putative roots up to this. Should be fmax/2.

%Define parameters 
Relationship_type_by_triplet_ID = obj.TripletTable.Relationship_Type(:);
freqs = obj.freqs;
R_max = obj.freqs(end)/2; %The maximum root frequencies considered is half the highest frequency.

%init data structures for results
obj.circR_analysis.R1_R2_Frequency_Indices = cell(obj.N_triplets,1);
obj.circR_analysis.standard_R1_R2 = NaN(length(0:0.5:R_max),length(0:0.5:R_max),20);
for triplet_ID = 1:obj.N_triplets
    Relationship_type = Relationship_type_by_triplet_ID(triplet_ID);
    XYZ_freqs = freqs(obj.circR_analysis.Frequency_Indices{Relationship_type});    
    data_corrected = obj.circR_analysis.Test.circR{Relationship_type}-obj.circR_analysis.Surrogate.circR{Relationship_type};
    data_uncorrected = obj.circR_analysis.Test.circR{Relationship_type};
    data_surrogate = obj.circR_analysis.Surrogate.circR{Relationship_type};
    x_R1 = obj.TripletTable.XY_R1(triplet_ID,1);
    y_R1 = obj.TripletTable.XY_R1(triplet_ID,2);
    x_R2 = obj.TripletTable.XY_R2(triplet_ID,1);
    y_R2 = obj.TripletTable.XY_R2(triplet_ID,2);
    R1_R2_freqs = [x_R1*XYZ_freqs(:,1)+y_R1*XYZ_freqs(:,2),x_R2*XYZ_freqs(:,1)+y_R2*XYZ_freqs(:,2)];
           
    good_inds = R1_R2_freqs(:,1)<=R_max & R1_R2_freqs(:,2)<=R_max & R1_R2_freqs(:,1) > R1_R2_freqs(:,2);
    R1_R2_freqs = R1_R2_freqs(good_inds,:);
    obj.circR_analysis.R1_R2_Frequency_Indices{triplet_ID} = R1_R2_freqs; %Correction: freqs, not inds
    
    goodData_corrected = data_corrected(good_inds);
    goodData_uncorrected = data_uncorrected(good_inds);
    goodData_surrogate = data_surrogate(good_inds);
    
    %for corrected
    if length(goodData_corrected)>1
        obj.circR_analysis.standard_R1_R2_corrected(:,:,triplet_ID) = ...
            griddata(R1_R2_freqs(:,1),R1_R2_freqs(:,2),goodData_corrected,obj.standard_grid_R1,obj.standard_grid_R2,'cubic'); %z4 will have NaN's where there's no data
    else
        obj.circR_analysis.standard_R1_R2_corrected(:,:,triplet_ID) = NaN(size(obj.standard_grid_R1));
    end
    
    %for uncorrected
    if length(goodData_uncorrected)>1
        obj.circR_analysis.standard_R1_R2_uncorrected(:,:,triplet_ID) = ...
            griddata(R1_R2_freqs(:,1),R1_R2_freqs(:,2),goodData_uncorrected,obj.standard_grid_R1,obj.standard_grid_R2,'cubic'); %z4 will have NaN's where there's no data
    else
        obj.circR_analysis.standard_R1_R2_uncorrected(:,:,triplet_ID) = NaN(size(obj.standard_grid_R1));
    end
    
    %for surrogate
    if length(goodData_surrogate)>1
        obj.circR_analysis.standard_R1_R2_surrogate(:,:,triplet_ID) = ...
            griddata(R1_R2_freqs(:,1),R1_R2_freqs(:,2),goodData_surrogate,obj.standard_grid_R1,obj.standard_grid_R2,'cubic'); %z4 will have NaN's where there's no data
    else
        obj.circR_analysis.standard_R1_R2_surrogate(:,:,triplet_ID) = NaN(size(obj.standard_grid_R1));
    end
        
end
end
