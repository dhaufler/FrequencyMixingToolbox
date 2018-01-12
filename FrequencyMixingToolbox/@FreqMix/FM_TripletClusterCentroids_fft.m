function [obj] = FM_TripletClusterCentroids_fft(obj,threshold,N_max)
%Use with output from 'WrappedPhaseTriplets.m'
%   Note:  X,Y,Z are indices, not frequencies


%init parameters
n_freqs     = obj.N_freqs;
freq_inds = obj.dF:0.25:n_freqs;

%output fields
obj.fft_pII_analysis.pII_XYZ_peaks = cell(obj.N_distinct_XYZ_relationships,1);

%ADOPT the following:
% Additional analysis, October 10th: Recover XYZ peaks.
XYZ_R = cell(8,1);
for freq_rel_type = 1:8
    v = obj.fft_pII_analysis.Test.pIIs{freq_rel_type}(:,1);
    X = obj.Frequency_Indices{freq_rel_type};
    F = scatteredInterpolant(X(:,1:2),v,'linear','none');
    %[x,y] = meshgrid(0.25:0.25:320,0.25:0.25:320);
    [x,y] = meshgrid(obj.dF:0.25:n_freqs,obj.dF:0.25:n_freqs); %create a mesh of interpolated INDICES
    v2___ = F([x(:),y(:)]);
    v2__ = reshape(v2___,size(x));
    v2_ = imgaussfilt(v2__,2);
    min_v2 = min(v2_(:));
    v2 = v2_-min_v2;
    [peaks_XY]=FastPeakFind(v2, threshold*2^16);
    peaks_X_ = freq_inds(peaks_XY(1:2:end));
    peaks_Y_ = freq_inds(peaks_XY(2:2:end));
    peaks_Z_ = [obj.distinct_XYZ_relationships(freq_rel_type,1)*peaks_X_ + obj.distinct_XYZ_relationships(freq_rel_type,2)*peaks_Y_]/(-obj.distinct_XYZ_relationships(freq_rel_type,3));
    IND = sub2ind(size(v2),peaks_XY(2:2:end),peaks_XY(1:2:end));
    %peaks_R_ = v2(IND)+min_v2;
    peaks_R_ = v2__(IND)-nanmedian(v2__(:));
    
    %sort, remove those that are too near to one another, and take up to 8
    [~,peaks_IND] = sort(peaks_R_,'descend');
    if length(peaks_IND)>0
    good_Inds_ = peaks_IND(1);
    n_good_inds = 1;
    for w=2:length(peaks_IND)
        if norm([peaks_X_(peaks_IND(w)) - peaks_X_(good_Inds_(n_good_inds)),...
                peaks_Y_(peaks_IND(w)) - peaks_Y_(good_Inds_(n_good_inds)),...
                peaks_Z_(peaks_IND(w)) - peaks_Z_(good_Inds_(n_good_inds))]) > 5
            good_Inds_ = [good_Inds_;peaks_IND(w)];
            n_good_inds = n_good_inds+1;
        end
    end

    good_Inds = good_Inds_(1:min([length(good_Inds_),N_max]));
    peaks_X = peaks_X_(good_Inds);
    peaks_Y = peaks_Y_(good_Inds);
    peaks_Z = peaks_Z_(good_Inds);
    peaks_R = peaks_R_(good_Inds);        
    XYZ_R{freq_rel_type,1} = [peaks_X',peaks_Y',peaks_Z',peaks_R];    
    else
    XYZ_R{freq_rel_type,1}=[];
    end
    end
obj.fft_pII_analysis.pII_XYZ_peaks  = XYZ_R; %changed this on Linux version


end

