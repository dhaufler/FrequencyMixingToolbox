function [obj] = image_S_pyramid(obj,S_threshold)

%init
[X_inds_,Y_inds_,Z_inds_] = ndgrid(obj.freqs,obj.freqs,obj.freqs);
X_inds = X_inds_+(rand(size(X_inds_))-0.5)*obj.dF;
Y_inds = Y_inds_+(rand(size(X_inds_))-0.5)*obj.dF;
Z_inds = Z_inds_+(rand(size(X_inds_))-0.5)*obj.dF;
                    
cluster_colors_8 = zeros(8,3);
cluster_colors_8(:,:) =   [68, 91, 232;...
                                     229, 43, 198;...
                                     232, 68, 68;...
                                     232, 134, 68;...
                                     232, 197,68;...
                                     205, 232, 68;...
                                     127, 232, 68;...
                                     68, 224, 232]/255;


h=figure; subplot(5,5,[7 8 9 12 13 14 17 18 19]); hold on;

%----this was the original-----
%pairwise_control = sum(obj.fft_pII_analysis.S_full(:,:,:,9:11),4);
pairwise_control = nthroot(obj.fft_pII_analysis.S_full(:,:,:,9),3).*nthroot(obj.fft_pII_analysis.S_full(:,:,:,10),3).*nthroot(obj.fft_pII_analysis.S_full(:,:,:,11),3);

% %----this is to set a threshold by percentile----
% S_all_clu_types = zeros(size(obj.fft_pII_analysis.S_full(:,:,:,1:8)));
% for i=1:8
%     S_all_clu_types(:,:,:,i) = obj.fft_pII_analysis.S_full(:,:,:,i)  - pairwise_control;
% end
% S_max_clu_type= max(S_all_clu_types,[],4); % consider only the largest spectral-pII value for any triplet
% temp = S_max_clu_type(1:120,1:120,1:120);
% S_threshold = prctile(temp(:),99)
% %S_threshold =0.2073;
for i=1:8
    S = obj.fft_pII_analysis.S_full(:,:,:,i)  - pairwise_control;
    %[~,max_ind] = max(obj.fft_pII_analysis.S_full,[],4);
    %i_is_max = max_ind==i;
    %plot_cluster(X_inds(i_is_max),Y_inds(i_is_max),Z_inds(i_is_max),S(i_is_max),squeeze(cluster_colors_8(i,:)),S_threshold,h)
    plot_cluster(X_inds,Y_inds,Z_inds,S,squeeze(cluster_colors_8(i,:)),S_threshold,h)
end
format_figure(h)
end
%-------------------------------

% %----try computing individual thresholds, and filtering-----
% pairwise_control = sum(obj.fft_pII_analysis.S_full(:,:,:,9:11),4);
% for i=1:8
%     S_ = obj.fft_pII_analysis.S_full(:,:,:,i)  - pairwise_control;
%     S_threshold = prctile(S_(:),98);
%     S= medfilt3(S_,[5 5 5]);
%     plot_cluster(X_inds,Y_inds,Z_inds,S,squeeze(cluster_colors_8(i,:)),S_threshold,h)
% end
% format_figure(h)
% end
% %-------------------------------


function h = plot_cluster(X_inds_,Y_inds_,Z_inds_,S_,c,S_threshold,h)
figure(h); subplot(5,5,[7 8 9 12 13 14 17 18 19]);
goodInds =  S_>S_threshold; 
S = S_(goodInds);
X_inds = X_inds_(goodInds);
Y_inds = Y_inds_(goodInds);
Z_inds = Z_inds_(goodInds);
y = quantile(S,[.33 .66]);

qt1 = S<y(1);
qt2 = S>=y(1)&S<y(2);
qt3 = S>y(2);


%plot in 4 quartiles
scatter3(X_inds(qt1),Y_inds(qt1),Z_inds(qt1),5 ,0.4*squeeze(c(:))','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0)
scatter3(X_inds(qt2),Y_inds(qt2),Z_inds(qt2),5 ,0.6*squeeze(c(:))','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0)
scatter3(X_inds(qt3),Y_inds(qt3),Z_inds(qt3),5 ,squeeze(c(:))','filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0)
end

function h = format_figure(h)
figure(h); subplot(5,5,[7 8 9 12 13 14 17 18 19]);
xlim([0 300])
ylim([0 200])
zlim([0 160])
xlabel('X'); ylabel('Y');zlabel('Z')
daspect([0.5,0.5,0.5]);
box on; grid on;
set(gca,'cameraviewanglemode','manual')
end