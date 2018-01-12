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
Relationship_type_by_triplet_ID = [1 1 1 1 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
freqs = obj.freqs;
R_max = obj.freqs(end)/2; %The maximum root frequencies considered is half the highest frequency.

%init data structures for results
obj.circR_analysis.R1_R2_Frequency_Indices = cell(20,1);
obj.circR_analysis.standard_R1_R2 = NaN(length(0:0.25:R_max),length(0:0.25:R_max),20);

for triplet_ID = 1:20
    Relationship_type = Relationship_type_by_triplet_ID(triplet_ID);
    XYZ_freqs = freqs(obj.circR_analysis.Frequency_Indices{Relationship_type});
    
    %data = obj.circR_analysis.Test.circR{Relationship_type};
    temp1 = obj.circR_analysis.Control.circR{Relationship_type};
    temp2 = obj.circR_analysis.Control.circRp{Relationship_type};
    
    %temp1(temp2>0.01)=0;
    pairwise_relationship = max(temp1); 
    data = obj.circR_analysis.Test.circR{Relationship_type};
    %data = obj.circR_analysis.Test.circR{Relationship_type}-pairwise_relationship;  
    
    switch triplet_ID
        case 1
            R1_R2_freqs = [XYZ_freqs(:,1),XYZ_freqs(:,3)]; 
        case 2
            R1_R2_freqs = [XYZ_freqs(:,1),XYZ_freqs(:,2)]; 
        case 3
            R1_R2_freqs = [XYZ_freqs(:,2),XYZ_freqs(:,3)]; 
        case 4
            R1_R2_freqs = [XYZ_freqs(:,1)-XYZ_freqs(:,2)/2,XYZ_freqs(:,2)/2]; 
        case 5
            R1_R2_freqs = [XYZ_freqs(:,1)-XYZ_freqs(:,3)/2,XYZ_freqs(:,3)/2]; 
        case 6
            R1_R2_freqs = [XYZ_freqs(:,1)/2,XYZ_freqs(:,2)-XYZ_freqs(:,1)/2]; 
        case 7
            R1_R2_freqs = [XYZ_freqs(:,1)-XYZ_freqs(:,3),XYZ_freqs(:,3)];        
        case 8
            R1_R2_freqs = [XYZ_freqs(:,1)/2,XYZ_freqs(:,2)/2];
        case 9
            R1_R2_freqs = [XYZ_freqs(:,1)-XYZ_freqs(:,2),XYZ_freqs(:,2)];
        case 10
            R1_R2_freqs = [XYZ_freqs(:,1)/2,XYZ_freqs(:,3)/2];
        case 11
            R1_R2_freqs = [XYZ_freqs(:,2),XYZ_freqs(:,1)-XYZ_freqs(:,2)];
        case 12
            R1_R2_freqs = [XYZ_freqs(:,1)/2,XYZ_freqs(:,3)/2];
        case 13
            R1_R2_freqs = [XYZ_freqs(:,1),XYZ_freqs(:,1)-XYZ_freqs(:,3)];
        case 14
            R1_R2_freqs = [XYZ_freqs(:,3),XYZ_freqs(:,2)/2];
        case 15
            R1_R2_freqs = [XYZ_freqs(:,1),XYZ_freqs(:,3)/2];
        case 16
            R1_R2_freqs = [XYZ_freqs(:,2),XYZ_freqs(:,3)/2];
        case 17
            R1_R2_freqs = [XYZ_freqs(:,1)/2,XYZ_freqs(:,2)];
        case 18
            R1_R2_freqs = [XYZ_freqs(:,1)/2,XYZ_freqs(:,3)];
        case 19
            R1_R2_freqs = [XYZ_freqs(:,2),XYZ_freqs(:,1)/2];
        case 20
            R1_R2_freqs = [XYZ_freqs(:,1)/2,XYZ_freqs(:,3)];
    end
    obj.circR_analysis.R1_R2_Frequency_Indices{triplet_ID} = R1_R2_freqs;

    %limit output to frequencies below R_max    
    good_inds = R1_R2_freqs(:,1)<=R_max & R1_R2_freqs(:,2)<=R_max;
    R1_R2_freqs = R1_R2_freqs(good_inds,:);
    goodData = data(good_inds);
    
    %interpolate result onto standard_grid
    obj.circR_analysis.standard_R1_R2(:,:,triplet_ID) = ...
        griddata(R1_R2_freqs(:,1),R1_R2_freqs(:,2),goodData,obj.standard_grid_R1,obj.standard_grid_R2,'cubic'); %z4 will have NaN's where there's no data

    
end




end
