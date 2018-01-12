function [boot_p] = FM_bootstrap1(R1,R2,rootsScore,X,freqs,Triplet_Table)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

%init
boot_p = ones(length(R1),1);
N_root_pairs = length(R1);
N_boot_samples = 100;
N_samples = size(X,1);
N_triplets = size(Triplet_Table,1); %should be 20
Component_Index =  [1 0; ... R1
                    0 1; ... R2
                    2 0; ... H1
                    0 2; ... H2
                    1 1; ... E1
                    1 -1];%   E2

X_ind = zeros(20,1); X_coef = zeros(20,1);
Y_ind = zeros(20,1); Y_coef = zeros(20,1);
Z_ind = zeros(20,1); Z_coef = zeros(20,1);
for i=1:N_triplets
    %better if this wasn't recomputed each time
    X_ind(i) = find(Component_Index(:,1)==Triplet_Table.compX(i,1) & Component_Index(:,2)==Triplet_Table.compX(i,2));
    Y_ind(i) = find(Component_Index(:,1)==Triplet_Table.compY(i,1) & Component_Index(:,2)==Triplet_Table.compY(i,2));
    Z_ind(i) = find(Component_Index(:,1)==Triplet_Table.compZ(i,1) & Component_Index(:,2)==Triplet_Table.compZ(i,2));
        
    X_coef(i) = Triplet_Table.XYZ_relationship(i,1);
    Y_coef(i) = Triplet_Table.XYZ_relationship(i,2);
    Z_coef(i) = Triplet_Table.XYZ_relationship(i,3);
end

for i=1:N_root_pairs
    bootstat = [];
    [FM_components,FM_freqs] = getComponents(R1(i),R2(i),X,freqs);
    S__ = FM_surrogate(FM_components./abs(FM_components),N_boot_samples*N_samples);
    S_ = angle(S__);
    for j=1:N_boot_samples
        
        %S_ = FM_surrogate(FM_components./abs(FM_components),size(FM_components,1));
        S = S_((j-1)*N_samples+1:j*N_samples,:);
        bootstat_temp = FM_Assess_Putative_Roots(S,FM_freqs,X_ind,Y_ind,Z_ind,X_coef,Y_coef,Z_coef);
        bootstat = [bootstat,bootstat_temp];
    end
    bootScore = nanmedian(bootstat,1);
    boot_p(i) = sum(bootScore>=rootsScore(i))/N_boot_samples;
end
end


