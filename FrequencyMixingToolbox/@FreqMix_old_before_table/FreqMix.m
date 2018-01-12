classdef FreqMix < dynamicprops
    %
    
    properties
        data        %complex matrix of samples-by-frequencies
        S_data      %surrogate data
        freqs       %vector of frequencies
        freqInds    %vector of frequency indices
        N_freqs     %number of frequencies
        N_samples   %number of samples
        dF          %frequency increment
        standard_grid_R1
        standard_grid_R2
        
        %below are structures storing analysis-specific parameters and
        %results
        pII_analysis    
        circR_analysis
        surrogate_analysis
    end
    methods
        
        % constructor/initialization
        function obj = FreqMix(data_In,freqs_In)
        %   data_In:    n-by-p complex numbers, where n is the number
        %               of observations, and p is the number of frequencies.
        %   freqs_In:   vector of p frequencies.
        
            %check that input dimensions are consistent
            if size(data_In,2)~=length(freqs_In)
                error('Input dimensions are inconsistent. Spectral input should be n-by-p, where n is the number of samples and p is the number of frequencies');
            end

            %init data and parameters
            obj.data = data_In;
            obj.S_data = FM_surrogate(data_In);
            obj.freqs = freqs_In;
            obj.N_freqs = size(data_In,2);
            obj.N_samples = size(data_In,1);
            obj.freqInds = 1:length(freqs_In);
            obj.dF = freqs_In(2) - freqs_In(1);
            
            %for interpolating R1-R2 estimates onto a standard grid
            [obj.standard_grid_R1, obj.standard_grid_R2] = meshgrid(0:0.25:(obj.freqs(end)/2));
            
            obj.circR_analysis = [];
            obj.FM_circR
            obj.FM_TripletClusterCentroids(0.2)
            obj.FM_R1_R2
        end
        % additional methods
        obj = FM_circR(obj)
        obj = FM_TripletClusterCentroids(obj,threshold)
        obj = FM_R1_R2(obj)
        
        %imaging and visualization
        obj = image_PhaseTriplet(obj,freq1,freq2,freq3) %plots a phase triplet distribution
        obj = image_XYZ_plane(obj,relationship_type) %plots relationship-type 1-8 in 3D frequency space
        obj = image_R1_R2_by_triplets(obj,triplet_number)
    end
end
        
            
        