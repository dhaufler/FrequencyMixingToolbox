function obj = FM_get_cluster_orientation(obj)
%UNTITLED6 Summary of this function goes here
%   Run this after FM_TripletClusterCentroids_fft

obj.fft_pII_analysis.pII_XYZ_peaks_eigs  = cell(obj.N_distinct_XYZ_relationships,1);  
obj.fft_pII_analysis.pII_XYZ_peaks_eigvec  = cell(obj.N_distinct_XYZ_relationships,1);  
for freq_rel_type = 1:8
    obj.fft_pII_analysis.pII_XYZ_peaks_eigs{freq_rel_type} = NaN(size(obj.fft_pII_analysis.pII_XYZ_peaks{freq_rel_type},1),3,2);
    obj.fft_pII_analysis.pII_XYZ_peaks_eigvec{freq_rel_type} = NaN(size(obj.fft_pII_analysis.pII_XYZ_peaks{freq_rel_type},1),3,2,2);
    %obj.fft_pII_analysis.pII_XYZ_R{freq_rel_type} = NaN(size(obj.fft_pII_analysis.pII_XYZ_peaks{freq_rel_type},1),3);
    %obj.fft_pII_analysis.pII_XYZ_p{freq_rel_type} = NaN(size(obj.fft_pII_analysis.pII_XYZ_peaks{freq_rel_type},1),3);
    for peak_ind = 1:size(obj.fft_pII_analysis.pII_XYZ_peaks{freq_rel_type},1)
        peak_xyz = obj.fft_pII_analysis.pII_XYZ_peaks{freq_rel_type}(peak_ind,1:3);
        
        peak_xyz_surround_ = NaN(9,9,9);

        for i=-4:4
            try
            S_i = obj.data(:,round(peak_xyz(1))+3*i);
            catch
            end
            for j=-4:4
                try                    
                S_j = obj.data(:,round(peak_xyz(2))+3*j);
                catch
                end
                for k = -4:4
                    try
                    S_k = obj.data(:,round(peak_xyz(3))+3*k); 
                    catch
                    end
                    try
                    H = triplet_FFT_pII([S_i,S_j,S_k]); %compute pII at the specified points                                       
                    %peak_xyz_surround_(i+5,j+5,k+5) = (H.pII).*(H.S_rel(freq_rel_type)./sum(H.S_rel(:))); %scale pII by the color index of the relationship type
                     peak_xyz_surround_(i+5,j+5,k+5) = (H.pII); %scale pII by the color index of the relationship type
                    catch
                    end
                end
            end
        end
                
        peak_xyz_surround_(isnan(peak_xyz_surround_)) = nanmean(peak_xyz_surround_(:));  %set any NaNs to mean      
        peak_xyz_surround = imgaussfilt3(peak_xyz_surround_,[1 1 1],'padding','replicate');
        %get subscripts of largest element in the computed volume
        [~,max_ind] = max(peak_xyz_surround(:));
        [Imax,Jmax,Kmax] = ind2sub(size(peak_xyz_surround),max_ind);

        % Compute the Hessian matrices along the different intersecting planes
        [gy,gx,gz] = gradient(peak_xyz_surround);
        [gxy,gxx,gxz] = gradient(gx);
        [gyy,gyx,gyz] = gradient(gy);
        [gzy,gzx,gzz] = gradient(gz);
        H_xy =  [gxx(Imax,Jmax,Kmax), gxy(Imax,Jmax,Kmax);...
                gyx(Imax,Jmax,Kmax), gyy(Imax,Jmax,Kmax)];
        H_xz =  [gxx(Imax,Jmax,Kmax), gxz(Imax,Jmax,Kmax);...
                gzx(Imax,Jmax,Kmax), gzz(Imax,Jmax,Kmax)];
        H_yz =  [gyy(Imax,Jmax,Kmax), gyz(Imax,Jmax,Kmax);...
                gzy(Imax,Jmax,Kmax), gzz(Imax,Jmax,Kmax)];                
                     
         [V_xy,D_xy] = eig(H_xy);
         [V_xz,D_xz] = eig(H_xz); 
         [V_yz,D_yz] = eig(H_yz);
         
         %just save V and D
         obj.fft_pII_analysis.pII_XYZ_peaks_eigs{freq_rel_type}(peak_ind,1,1:2) = [D_xy(1,1),D_xy(2,2)];
         obj.fft_pII_analysis.pII_XYZ_peaks_eigs{freq_rel_type}(peak_ind,2,1:2) = [D_xz(1,1),D_xz(2,2)];
         obj.fft_pII_analysis.pII_XYZ_peaks_eigs{freq_rel_type}(peak_ind,3,1:2) = [D_yz(1,1),D_yz(2,2)];
         
         obj.fft_pII_analysis.pII_XYZ_peaks_eigvec{freq_rel_type}(peak_ind,1,:,:) = V_xy;
         obj.fft_pII_analysis.pII_XYZ_peaks_eigvec{freq_rel_type}(peak_ind,2,:,:) = V_xz;
         obj.fft_pII_analysis.pII_XYZ_peaks_eigvec{freq_rel_type}(peak_ind,3,:,:) = V_yz;
        
    end    
end
        

end

