function [xy_corr, xz_corr, yz_corr] = xyz_ellipsoid_corr(X)

X_filt = imgaussfilt3(X,[2 2 2],'padding','replicate');
[~,max_ind] = max(X_filt(:));
[Imax,Jmax,Kmax] = ind2sub(size(X_filt),max_ind);

[gy,gx,gz] = gradient(X_filt);
[gxy,gxx,gxz] = gradient(gx);
[gyy,gyx,gyz] = gradient(gy);
[gzy,gzx,gzz] = gradient(gz);

H_xy = [ gxx(Imax,Jmax,Kmax), gxy(Imax,Jmax,Kmax);...
        gyx(Imax,Jmax,Kmax),  gyy(Imax,Jmax,Kmax)];

H_xz = [gxx(Imax,Jmax,Kmax), gxz(Imax,Jmax,Kmax);...
        gzx(Imax,Jmax,Kmax), gzz(Imax,Jmax,Kmax)];
    
H_yz = [gyy(Imax,Jmax,Kmax), gyz(Imax,Jmax,Kmax);... 
    	gzy(Imax,Jmax,Kmax), gzz(Imax,Jmax,Kmax)];

[V_xy,D_xy] = eig(H_xy);
[V_xz,D_xz] = eig(H_xz);
[V_yz,D_yz] = eig(H_yz);    

xy_corr_1 = V_xy(:,1)/D_xy(1,1);
xy_corr_2 = V_xy(:,2)/D_xy(2,2);

xz_corr_1 = V_xz(:,1)/D_xz(1,1);
xz_corr_2 = V_xz(:,2)/D_xz(2,2);

yz_corr_1 = V_yz(:,1)/D_yz(1,1);
yz_corr_2 = V_yz(:,2)/D_yz(2,2);

xy_corr = xy_corr_2/(norm(xy_corr_1));
xz_corr = xz_corr_2/(norm(xz_corr_1));
yz_corr = yz_corr_2/(norm(yz_corr_1));
end