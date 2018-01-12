function  [boot_p] = FM_bootstrap2(R1,R2,rootsScore,X,freqs,Triplet_Table,boot_p1)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
%init
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
    if boot_p1(i)<=0.01
        disp([num2str(i),' of ',num2str(N_root_pairs)])
        bootstat = [];
        [FM_components,FM_freqs] = getComponents(R1(i),R2(i),X,freqs);
        S__ = getFMcomponent_surrogate(FM_components./abs(FM_components));
        S_ = angle(S__);
        for j=1:N_boot_samples        
            S = S_(randsample(size(S_,1),N_samples,1),:);
            bootstat_temp = FM_Assess_Putative_Roots(S,FM_freqs,X_ind,Y_ind,Z_ind,X_coef,Y_coef,Z_coef);
            bootstat = [bootstat,bootstat_temp];
        end
        bootScore = nanmedian(bootstat,1);
        boot_p(i) = sum(bootScore>=rootsScore(i))/N_boot_samples;
    else
        boot_p(i) = NaN;
    end
end
end


function S = getFMcomponent_surrogate(Data_Matrix)
h=[];
d_scale = 2;
Max_iterations = 500;
%N_samples = 1000;

S = [];
for j=1:10
    S_start_ = FM_surrogate(Data_Matrix./abs(Data_Matrix),500);
    S_start = S_start_./abs(S_start_);
    new_S = getFMsurrogate3(Data_Matrix./abs(Data_Matrix),S_start,[1 1 0 0 0 0],[1 2], 0,1000,d_scale,h);
    for i=1:6
        new_S = getFMsurrogate3(Data_Matrix./abs(Data_Matrix),new_S, [1 1 1 1 1], i, 0, Max_iterations,d_scale,h);
    end
    S = [S;new_S];
end
end
