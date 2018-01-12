function [] = imageWPT(WPT_data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

f=WPT_data.f;
n_freqs = length(f);
n_tests = size(WPT_data.xyz_key,1);
p_c = 0.00001;


%change this section: I know which relationships correspond to which colors
color_index = zeros(n_tests,1);
for i=1:n_tests
    if ismember(sum(abs(WPT_data.xyz_key(i,:))),[2 3 4 5]) & min(abs(WPT_data.xyz_key(i,:)))==0
        color_index(i,1) = 2;
    else 
        color_index(i,1) = sum(abs(WPT_data.xyz_key(i,:)));
    end
end

%output values will only be set for X>Y>Z
good_vox = zeros(n_freqs,n_freqs,n_freqs);
for i=1:n_freqs
    for j=1:i-1
        for k=1:j-1
            good_vox(i,j,k) =1;
        end
    end
end

%select good blue-group
B_group_CRp_ = WPT_data.phaseCRp(:,:,:,color_index==3);
R_group_CRp_ = WPT_data.phaseCRp(:,:,:,color_index==4);
G_group_CRp_ = WPT_data.phaseCRp(:,:,:,color_index==5);
Y_group_CRp_ = WPT_data.phaseCRp(:,:,:,color_index==2);

B_group_CRp = min(B_group_CRp_,[],4);
R_group_CRp = min(R_group_CRp_,[],4);
G_group_CRp = min(G_group_CRp_,[],4);
Y_group_CRp = min(Y_group_CRp_,[],4);

%plot significant blue but not significant yellow
figure; hold on

B_group_indices = find(B_group_CRp<p_c & good_vox==1 & Y_group_CRp>=p_c);
R_group_indices = find(R_group_CRp<p_c & good_vox==1 & Y_group_CRp>=p_c);
G_group_indices = find(G_group_CRp<p_c & good_vox==1 & Y_group_CRp>=p_c);

% %like above, but without excluding yellow
% B_group_indices = find(B_group_CRp<p_c & good_vox==1);
% R_group_indices = find(R_group_CRp<p_c & good_vox==1);
% G_group_indices = find(G_group_CRp<p_c & good_vox==1);

Y_group_indices = find(Y_group_CRp<p_c & good_vox==1 & [(B_group_CRp<p_c)|(R_group_CRp<p_c)|(G_group_CRp<p_c)]);

[B1, B2, B3] = ind2sub(size(squeeze(WPT_data.phaseCRp(:,:,:,1))),B_group_indices);
[R1, R2, R3] = ind2sub(size(squeeze(WPT_data.phaseCRp(:,:,:,1))),R_group_indices);
[G1, G2, G3] = ind2sub(size(squeeze(WPT_data.phaseCRp(:,:,:,1))),G_group_indices);
[Y1, Y2, Y3] = ind2sub(size(squeeze(WPT_data.phaseCRp(:,:,:,1))),Y_group_indices);

B_Hx = f(B1)+2*(rand(size(B1'))-0.5);
B_Hy = f(B2)+2*(rand(size(B1'))-0.5);
B_Hz = f(B3)+2*(rand(size(B1'))-0.5);

R_Hx = f(R1)+2*(rand(size(R1'))-0.5);
R_Hy = f(R2)+2*(rand(size(R1'))-0.5);
R_Hz = f(R3)+2*(rand(size(R1'))-0.5);

G_Hx = f(G1)+2*(rand(size(G1'))-0.5);
G_Hy = f(G2)+2*(rand(size(G1'))-0.5);
G_Hz = f(G3)+2*(rand(size(G1'))-0.5);

Y_Hx = f(Y1)+2*(rand(size(Y1'))-0.5);
Y_Hy = f(Y2)+2*(rand(size(Y1'))-0.5);
Y_Hz = f(Y3)+2*(rand(size(Y1'))-0.5);

scatter3(B_Hx,B_Hy,B_Hz,30,'b','filled','MarkerFaceAlpha',0.3);%,10,all_comps_color,'filled','MarkerFaceAlpha',0.35, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.0);
scatter3(R_Hx,R_Hy,R_Hz,30,'r','filled','MarkerFaceAlpha',0.3);%,10,all_comps_color,'filled','MarkerFaceAlpha',0.35, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.0);
scatter3(G_Hx,G_Hy,G_Hz,30,'g','filled','MarkerFaceAlpha',0.3);%,10,all_comps_color,'filled','MarkerFaceAlpha',0.35, 'MarkerEdgeColor','k','MarkerEdgeAlpha',0.0);
%scatter3(Y_Hx,Y_Hy,Y_Hz,30,'y','filled','MarkerFaceAlpha',0.3);

daspect([0.5,0.5,0.5]);
set(gca,'cameraviewanglemode','manual')
box on; axis square; grid on;
xlim([f(1) f(end)]);
ylim([f(1) f(end)]);
zlim([f(1) f(end)]);
xlabel('X')
ylabel('Y')
zlabel('Z')
drawnow

end

