function [circR_dat] = FM_Assess_Putative_Roots(Component_Data,Component_Freqs,X_ind,Y_ind,Z_ind,X_coef,Y_coef,Z_coef)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

%check the input dimensions
if size(Component_Data,2)~=length(Component_Freqs)
    error('The number of data columns does not equal the number of specified frequencies')
end

%check component frequencies are correct
Component_Freqs = round(Component_Freqs);
if  Component_Freqs(3)    >   Component_Freqs(1)            ||...
    Component_Freqs(2)    ~=  2*Component_Freqs(1)          ||...
    Component_Freqs(4)    ~=  2*Component_Freqs(3)          ||...
    Component_Freqs(5)    ~=  Component_Freqs(1)+Component_Freqs(3)   ||...
    Component_Freqs(6)    ~=  Component_Freqs(1)-Component_Freqs(3)
    error('Component frequencies inconsistent. Order: R1,H1,R2,H2,E1,E2')
end


%init
N_triplets = 20;%size(Triplet_Table,1);
circR_dat = zeros(N_triplets,1);

% %frequencies of putative roots
% R1_f = Component_Freqs(1);
% R2_f = Component_Freqs(3); 

% %spectral data: complex-valued column vectors
% R1_dat = data(:,1); %data for frequency R1_f
% H1_dat = data(:,2); %data for frequency 2*R1_f
% R2_dat = data(:,3); %data for frequency R2_f
% H2_dat = data(:,4); %data for frequency 2*R2_f
% E1_dat = data(:,5); %data for frequency R1_f + R2_f 
% E2_dat = data(:,6); %data for frequency R1_f - R2_f

% For indexing by the coefficients of (R1,R2) in Triplet_Table
% Component_Index =  [1 0; ... R1
%                     0 1; ... R2
%                     2 0; ... H1
%                     0 2; ... H2
%                     1 1; ... E1
%                     1 -1];%   E2
%Component_Data =    [R1_dat R2_dat H1_dat H2_dat E1_dat E2_dat];

% Component_Freqs =   [R1_f,      ...
%                      R2_f,      ...
%                      2*R1_f,    ...
%                      2*R2_f,    ...
%                      R1_f+R2_f, ...
%                      R1_f-R2_f];

%determine which of the 20 triplets apply
%triplet_mask = zeros(20,1);
for i=1:20
%     %better if this wasn't recomputed each time
%     X_ind = find(Component_Index(:,1)==Triplet_Table.compX(i,1) & Component_Index(:,2)==Triplet_Table.compX(i,2));
%     Y_ind = find(Component_Index(:,1)==Triplet_Table.compY(i,1) & Component_Index(:,2)==Triplet_Table.compY(i,2));
%     Z_ind = find(Component_Index(:,1)==Triplet_Table.compZ(i,1) & Component_Index(:,2)==Triplet_Table.compZ(i,2));
    
    X_f = Component_Freqs(X_ind(i));
    Y_f = Component_Freqs(Y_ind(i));
    Z_f = Component_Freqs(Z_ind(i));
    
    if ((X_f > Y_f) && (X_f > Z_f) && (Y_f > Z_f))
        X_dat = Component_Data(:,X_ind(i));
        Y_dat = Component_Data(:,Y_ind(i));
        Z_dat = Component_Data(:,Z_ind(i));
        TestDistribution = mod((X_dat)*X_coef(i) + (Y_dat)*Y_coef(i) + (Z_dat)*Z_coef(i) + pi,2*pi)-pi;
        circR_dat(i) = circ_r(TestDistribution);
    else
        circR_dat(i,:) = NaN;
    end
end
               
%         X_coef = Triplet_Table.XYZ_relationship(i,1);
%         Y_coef = Triplet_Table.XYZ_relationship(i,2);
%         Z_coef = Triplet_Table.XYZ_relationship(i,3);
        
%         TestDistribution = mod(X_coef(i)*angle(X_dat) + Y_coef(i)*angle(Y_dat) + Z_coef(i)*angle(Z_dat) + pi,2*pi)-pi;
%         circR_dat(i,1) = circ_r(TestDistribution);
%         %circR_dat(i,2) = circ_rtest(TestDistribution);
% 
%     
% X_dat = Component_Data(:,X_ind);
% Y_dat = Component_Data(:,Y_ind);
% Z_dat = Component_Data(:,Z_ind);
% TestDistribution = mod(angle(X_dat)*repmat(X_coef,[1 20]) + angle(Y_dat)*repmat(Y_coef,[1 20]) + angle(Z_dat)*repmat(Z_coef,[1 20]) + pi,2*pi)-pi;
% circR_dat = circ_r(TestDistribution);
% circR_dat(logical(triplet_mask)) = NaN;
%circR_dat(i,2) = circ_rtest(TestDistribution);
end

