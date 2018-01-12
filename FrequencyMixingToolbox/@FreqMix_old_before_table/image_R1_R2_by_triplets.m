function obj = image_R1_R2_by_triplets(obj,triplet_IDs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

temp_data = obj.circR_analysis.standard_R1_R2(:,:,[triplet_IDs]);
%temp_data(temp_data<0)=0;

figure; 

if length(triplet_IDs)>1
    obj.circR_analysis.standard_Rvals_mean = nanmean(temp_data,3);
    obj.circR_analysis.standard_Rvals_median = nanmedian(temp_data,3);
    
    subplot(2,2,1)
    
    max_roots = max([obj.circR_analysis.standard_Rvals_mean(:); obj.circR_analysis.standard_Rvals_median(:)]);
    %colormap(bluewhitered)
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,obj.circR_analysis.standard_Rvals_mean); 
    set(h,'edgecolor','none'); caxis([-0.2 0.2]);colormap(bluewhitered)%caxis([0 max_roots]); 
    axis square; %colorbar
    freezeColors
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title('mean')

    subplot(2,2,3)
    %colormap(bluewhitered)
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,obj.circR_analysis.standard_Rvals_median); 
    set(h,'edgecolor','none'); caxis([-0.2 0.2]);%colormap(bluewhitered)%caxis([0 max_roots]); 
    axis square; %colorbar
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title('median')
    freezeColors
    colormap(parula)
    subplot(2,2,[2 4])
    imagesc(1:size(obj.data,1),obj.freqs,log(abs(obj.data'))); axis xy; ylim([0 obj.freqs(end)/2])
    ylabel('frequency (Hz)')
    xlabel('samples')
    suptitle(['Average of R1 vs R2 for average of triplets: ',num2str(triplet_IDs)])
    %suptitle(['For i =  ',num2str(image_number)])
else
    standard_Rvals =temp_data;
    max_roots = max(standard_Rvals(:));
    h = pcolor(obj.standard_grid_R1,obj.standard_grid_R2,standard_Rvals); 
    set(h,'edgecolor','none'); caxis([0 max_roots]); axis square; colorbar
    xlabel('R1 (Hz)');
    ylabel('R2 (Hz)');
    title(['R1 vs R2 for triplet: ',num2str(triplet_IDs)])
end

end

