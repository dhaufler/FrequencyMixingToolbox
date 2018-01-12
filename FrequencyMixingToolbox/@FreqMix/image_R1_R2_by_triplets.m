function obj = image_R1_R2_by_triplets(obj,triplet_IDs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% temp_data_corrected = obj.circR_analysis.standard_R1_R2_corrected(:,:,[triplet_IDs]);
% temp_data_uncorrected = obj.circR_analysis.standard_R1_R2_uncorrected(:,:,[triplet_IDs]);
% temp_data_surrogate = obj.circR_analysis.standard_R1_R2_surrogate(:,:,[triplet_IDs]);

temp_data_corrected     = obj.fft_pII_analysis.standard_R1_R2_corrected(:,:,[triplet_IDs]);
temp_data_uncorrected   = obj.fft_pII_analysis.standard_R1_R2_uncorrected(:,:,[triplet_IDs]);
temp_data_surrogate     = obj.fft_pII_analysis.standard_R1_R2_surrogate(:,:,[triplet_IDs]);

%temp_data_corrected(isnan(temp_data_corrected))=0;
%temp_data_uncorrected(isnan(temp_data_uncorrected))=0;
%temp_data_surrogate(isnan(temp_data_surrogate))=0;

figure; 

if length(triplet_IDs)>0
    %obj.circR_analysis.standard_Rvals_mean = nanmean(temp_data_corrected,3);
    %obj.circR_analysis.standard_Rvals_median = nanmedian(temp_data_corrected,3); 
    
    colormap(parula)
    subplot(2,6,[6 12])
    imagesc(1:size(obj.data,1),obj.freqs,log(abs(obj.data'))); axis xy; ylim([0 obj.freqs(end)/2]);
    freezeColors
    ylabel('frequency (Hz)')
    xlabel('samples')
    
    suptitle(['Average of R1 vs R2 for average of triplets: ',num2str(triplet_IDs)])
    %suptitle(['For i =  ',num2str(image_number)])
    
    subplot(2,6,1)    
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,nanmean(temp_data_corrected,3)); 
    set(h,'edgecolor','none'); caxis([-0.2 0.2]);colormap(bluewhitered)%caxis([0 max_roots]); 
    axis square; %colorbar
    freezeColors
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title('mean, corrected')
%     
    subplot(2,6,2) 
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,nanmean(temp_data_uncorrected,3)); 
    set(h,'edgecolor','none'); caxis([-0.2 0.2]);colormap(bluewhitered)%caxis([0 max_roots]); 
    axis square; %colorbar
    freezeColors
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title('mean, uncorrected')   

    subplot(2,6,3)
    %colormap(bluewhitered)
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,nanmean(temp_data_surrogate,3)); 
    set(h,'edgecolor','none'); caxis([-0.2 0.2]);%colormap(bluewhitered)%caxis([0 max_roots]); 
    axis square; %colorbar
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title('mean, surrogate')
    freezeColors
    
    subplot(2,6,7)    
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,nanmedian(temp_data_corrected,3)); 
    set(h,'edgecolor','none'); caxis([-0.01 0.01]);colormap(bluewhitered)%caxis([0 max_roots]); 
    axis square; %colorbar
    freezeColors
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title('median, corrected')
    
    subplot(2,6,8)    
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,nanmedian(temp_data_uncorrected,3)); 
    set(h,'edgecolor','none'); caxis([-0.01 0.01]);colormap(bluewhitered)%caxis([0 max_roots]); 
    axis square; %colorbar
    freezeColors
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title('median, uncorrected')

    subplot(2,6,9)    
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,nanmedian(temp_data_surrogate,3)); 
    set(h,'edgecolor','none'); caxis([-0.01 0.01]);colormap(bluewhitered)%caxis([0 max_roots]); 
    axis square; %colorbar
    freezeColors
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title('median, surrogate')

% else
%     standard_Rvals =temp_data_corrected;
%     max_roots = max(standard_Rvals(:));
%     h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,standard_Rvals); 
%     set(h,'edgecolor','none'); caxis([0 max_roots]); axis square; colorbar
%     xlabel('R1 (Hz)');
%     ylabel('R2 (Hz)');
%     title(['R1 vs R2 for triplet: ',num2str(triplet_IDs)])
% end

end

