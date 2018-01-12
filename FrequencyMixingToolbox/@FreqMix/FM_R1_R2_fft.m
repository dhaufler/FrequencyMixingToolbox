function [obj] = FM_R1_R2_fft(obj)

%Define parameters 
Relationship_type_by_triplet_ID = obj.TripletTable.Relationship_Type(:);
freqs = obj.freqs;
R_max = obj.freqs(end)/2; %The maximum root frequencies considered is half the highest frequency.

%init data structures for results
obj.fft_pII_analysis.R1_R2_Frequency_Indices = cell(obj.N_triplets,1);
obj.fft_pII_analysis.standard_R1_R2 = NaN(length(0:0.5:R_max),length(0:0.5:R_max),20);
for triplet_ID = 1:obj.N_triplets
    triplet_fft_inds = [(3+obj.TripletTable.XYZ_relationship(triplet_ID,1)),(3+ obj.TripletTable.XYZ_relationship(triplet_ID,2)),(3+ obj.TripletTable.XYZ_relationship(triplet_ID,3))];
    %I need the fftn subscripts for each of the relationships; these will
    %be derived from values in triplet table.
    Relationship_type = Relationship_type_by_triplet_ID(triplet_ID);
    XYZ_freqs = freqs(obj.Frequency_Indices{Relationship_type});
    
    data_uncorrected    = obj.fft_pII_analysis.Test.S{Relationship_type}(:,triplet_fft_inds(1),triplet_fft_inds(2),triplet_fft_inds(3));
    data_surrogate      = obj.fft_pII_analysis.Surrogate.S{Relationship_type}(:,triplet_fft_inds(1),triplet_fft_inds(2),triplet_fft_inds(3));
    data_corrected      = obj.fft_pII_analysis.Test.S{Relationship_type}(:,triplet_fft_inds(1),triplet_fft_inds(2),triplet_fft_inds(3))-...
                            obj.fft_pII_analysis.Surrogate.S{Relationship_type}(:,triplet_fft_inds(1),triplet_fft_inds(2),triplet_fft_inds(3));
    
    x_R1 = obj.TripletTable.XY_R1(triplet_ID,1);
    y_R1 = obj.TripletTable.XY_R1(triplet_ID,2);
    x_R2 = obj.TripletTable.XY_R2(triplet_ID,1);
    y_R2 = obj.TripletTable.XY_R2(triplet_ID,2);
    R1_R2_freqs = [x_R1*XYZ_freqs(:,1)+y_R1*XYZ_freqs(:,2),x_R2*XYZ_freqs(:,1)+y_R2*XYZ_freqs(:,2)];
           
    good_inds = R1_R2_freqs(:,1)<=R_max & R1_R2_freqs(:,2)<=R_max & R1_R2_freqs(:,1) > R1_R2_freqs(:,2);
    R1_R2_freqs = R1_R2_freqs(good_inds,:);
    obj.fft_pII_analysis.R1_R2_Frequency_Indices{triplet_ID} = R1_R2_freqs; %Correction: freqs, not inds
    
    goodData_uncorrected = data_uncorrected(good_inds);
    goodData_surrogate = data_surrogate(good_inds);
    goodData_corrected = data_corrected(good_inds);

    
    %for uncorrected
    if length(goodData_uncorrected)>1
        obj.fft_pII_analysis.standard_R1_R2_uncorrected(:,:,triplet_ID) = ...
            griddata(R1_R2_freqs(:,1),R1_R2_freqs(:,2),goodData_uncorrected,obj.standard_grid_R1,obj.standard_grid_R2,'cubic'); %z4 will have NaN's where there's no data
    else
        obj.fft_pII_analysis.standard_R1_R2_uncorrected(:,:,triplet_ID) = NaN(size(obj.standard_grid_R1));
    end
    
    %for surrogate
    if length(goodData_surrogate)>1
        obj.fft_pII_analysis.standard_R1_R2_surrogate(:,:,triplet_ID) = ...
            griddata(R1_R2_freqs(:,1),R1_R2_freqs(:,2),goodData_surrogate,obj.standard_grid_R1,obj.standard_grid_R2,'cubic'); %z4 will have NaN's where there's no data
    else
        obj.fft_pII_analysis.standard_R1_R2_surrogate(:,:,triplet_ID) = NaN(size(obj.standard_grid_R1));
    end
    
    %for corrected
    if length(goodData_corrected)>1
        obj.fft_pII_analysis.standard_R1_R2_corrected(:,:,triplet_ID) = ...
            griddata(R1_R2_freqs(:,1),R1_R2_freqs(:,2),goodData_corrected,obj.standard_grid_R1,obj.standard_grid_R2,'cubic'); %z4 will have NaN's where there's no data
    else
        obj.fft_pII_analysis.standard_R1_R2_corrected(:,:,triplet_ID) = NaN(size(obj.standard_grid_R1));
    end
%     
%     

        
end
end
