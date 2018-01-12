classdef FreqMix < dynamicprops
    %
    
    properties
        data        %complex matrix of samples-by-frequencies
        TripletTable %table of triplet types and relationships
        S_data      %surrogate data
        freqs       %vector of frequencies
        freqInds    %vector of frequency indices
        N_freqs     %number of frequencies
        N_samples   %number of samples
        dF          %frequency increment
        standard_grid_R1
        standard_grid_R2
        distinct_XYZ_relationships
        N_distinct_XYZ_relationships
        N_triplets
        Frequency_Indices
        
        %below are structures storing analysis-specific parameters and
        %results
        fft_pII_analysis        
        pII_analysis
        circR_analysis
        surrogate_analysis
    end
    methods
        
        % constructor/initialization
        function obj = FreqMix(data_In,freqs_In,Triplet_Table)
        %   data_In:    n-by-p complex numbers, where n is the number
        %               of observations, and p is the number of frequencies.
        %   freqs_In:   vector of p frequencies.
        
            %check that input dimensions are consistent
            if size(data_In,2)~=length(freqs_In)
                error('Input dimensions are inconsistent. Spectral input should be n-by-p, where n is the number of samples and p is the number of frequencies');
            end

            %init data and parameters
            obj.data = data_In;
            obj.TripletTable = Triplet_Table;
            obj.S_data = FM_surrogate(data_In,size(data_In,1));
            obj.freqs = freqs_In;
            obj.N_freqs = size(data_In,2);
            obj.N_samples = size(data_In,1);
            obj.freqInds = 1:length(freqs_In);
            obj.dF = freqs_In(2) - freqs_In(1);
            obj.distinct_XYZ_relationships = unique(Triplet_Table.XYZ_relationship(:,:),'rows','stable');
            obj.N_distinct_XYZ_relationships = size(obj.distinct_XYZ_relationships,1);
            obj.N_triplets = size(Triplet_Table,1);
            %for interpolating R1-R2 estimates onto a standard grid
            [obj.standard_grid_R1, obj.standard_grid_R2] = meshgrid(0:0.5:(obj.freqs(end)/2));
            obj.fft_pII_analysis = [];
            
            % These are the analyses run for Oct. 17 batch analysis
            obj.FM_init_XYZ_relationships;  
            obj.FM_FFT_pII;
            %obj.FM_TripletClusterCentroids_fft(0,8) %inputs are threshold, and maximum number of clusters
            %obj.FM_get_cluster_orientation
            %obj.circR_analysis = [];
            %obj.FM_circR
            %obj.FM_TripletClusterCentroids(0.2)
            %obj.FM_R1_R2
        end
        % additional methods
        obj = FM_init_XYZ_relationships(obj)
        obj = FM_FFT_pII(obj)
        obj = FM_R1_R2_fft(obj)
        obj = FM_TripletClusterCentroids_fft(obj,threshold,N_max)
        obj = FM_get_cluster_orientation(obj)
        obj = FM_FFT_S_full(obj)
                
        %old approach
        %obj = FM_circR(obj)
        %obj = FM_R1_R2(obj)
        %obj = FM_TripletClusterCentroids(obj,threshold)               
        %obj = FM_FFT_pII_complete(obj)          
        
        %imaging and visualization
        obj = image_S_pyramid(obj,S_threshold)
        obj = image_PhaseTriplet(obj,freq1,freq2,freq3) %plots a phase triplet distribution
        obj = image_XYZ_plane(obj,relationship_type) %plots relationship-type 1-8 in 3D frequency space
        obj = image_R1_R2_by_triplets(obj,triplet_number)
    end
end
        
            
        