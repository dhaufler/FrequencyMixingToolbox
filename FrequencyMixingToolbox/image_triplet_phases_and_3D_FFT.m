function [hh] = image_triplet_phases_and_3D_FFT(S)
% See 'Frequency_Mixing_Demo.m' to create a signal

%parameters and init
X = angle(S);
N_voxel_edges = 5; %use a power of 2 for FFT
N_phase_samples = size(S,1);
hh=figure;

%get binned histogram counts
N_XYZ_ = histcn(X, linspace(-pi,pi,N_voxel_edges+1), linspace(-pi,pi,N_voxel_edges+1), linspace(-pi,pi,N_voxel_edges+1));
N_XYZ = (N_XYZ_- mean(N_XYZ_))/sum(N_XYZ_(:));

% Compute 3D FFT
S_ = fftn(N_XYZ,size(N_XYZ)); S = (abs(S_)).^2;
[X_inds,Y_inds,Z_inds] = ndgrid(1:N_voxel_edges,1:N_voxel_edges,1:N_voxel_edges);   %indices of S

% 1. 3D scatterplot of phase angles
subplot(2,2,1)
scatter3(X(:,1),X(:,2),X(:,3),2,'k','filled'); 
xlim([-pi pi]);ylim([-pi pi]);zlim([-pi pi]);
xlabel('X'); ylabel('Y');zlabel('Z')
daspect([0.5,0.5,0.5]);
box on; axis square; grid on;
set(gca,'cameraviewanglemode','manual')

% 2. 3D histogram of binned phase angles
subplot(2,2,2); colormap(jet)
scatter3(X_inds(:),Y_inds(:),Z_inds(:),8000*[abs(N_XYZ(:))+0.00005],8000*[abs(N_XYZ(:))+0.00005],'filled','MarkerFaceAlpha',0.6)
xlim([1 5]); ylim([1 5]); zlim([1 5]);
xlabel('X'); ylabel('Y');zlabel('Z')
set(gca,'XTick',[1 2 3 4 5])
set(gca,'XTickLabel',{'-pi' '-pi/2' '0' 'pi/2' 'pi'});
set(gca,'YTick',[1 2 3 4 5])
set(gca,'YTickLabel',{'-pi' '-pi/2' '0' 'pi/2' 'pi'});
set(gca,'ZTick',[1 2 3 4 5])
set(gca,'ZTickLabel',{'-pi' '-pi/2' '0' 'pi/2' 'pi'});
daspect([0.5,0.5,0.5]);
box on; axis square; grid on;
set(gca,'cameraviewanglemode','manual')

% 3. 3D FFT, non-shifted
subplot(2,2,3); colormap(jet)
scatter3(X_inds(:),Y_inds(:),Z_inds(:),500*[S(:)+0.002],500*[S(:)+0.002],'filled','MarkerFaceAlpha',0.6)
xlim([1 5]); ylim([1 5]); zlim([1 5]);
xlabel('X'); ylabel('Y');zlabel('Z')
set(gca,'XTick',[1 2 3 4 5])
set(gca,'XTickLabel',[0 1 2 -2 -1]);
set(gca,'YTick',[1 2 3 4 5])
set(gca,'YTickLabel',[0 1 2 -2 -1]);
set(gca,'ZTick',[1 2 3 4 5])
set(gca,'ZTickLabel',[0 1 2 -2 -1]);
daspect([0.5,0.5,0.5]);
box on; axis square; grid on;
set(gca,'cameraviewanglemode','manual')

% 4. 3D FFT, shifted
S_shifted = fftshift(S);
subplot(2,2,4); colormap(jet)
scatter3(X_inds(:),Y_inds(:),Z_inds(:),500*[S_shifted(:)+0.00001],500*[S_shifted(:)+0.00001],'filled','MarkerFaceAlpha',0.6)
xlim([1 5]); ylim([1 5]); zlim([1 5]);
xlabel('X'); ylabel('Y');zlabel('Z')
set(gca,'XTick',[1 2 3 4 5])
set(gca,'XTickLabel',[-2 -1 0 1 2]);
set(gca,'YTick',[1 2 3 4 5])
set(gca,'YTickLabel',[-2 -1 0 1 2]);
set(gca,'ZTick',[1 2 3 4 5])
set(gca,'ZTickLabel',[-2 -1 0 1 2]);
daspect([0.5,0.5,0.5]);
box on; axis square; grid on;
set(gca,'cameraviewanglemode','manual')


% This doesn't seem to work when used within the function; execute at
% command window:
h1 = subplot(2,2,1);
h2 = subplot(2,2,2);
h3 = subplot(2,2,3);
h4 = subplot(2,2,4);
linkprop([h1 h2 h3 h4], 'View'); %h is the axes handle



end

