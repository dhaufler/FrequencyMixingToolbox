function [R1,R2,Ramp,pR1R2,D,D2,D3,D4,trp_type] = get_possible_R1_R2_from_set(X,Y,RelType,Centroid_mag, TripletTable)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

max_dist = 10;
minimum_agreeing_clusters = 4;
N_clusters = length(X);%total number of XYZ centroids for the microstate
pR1R2 = cell(N_clusters,1); %cell array of possible roots for the N clusters

for j=1:length(X) %for each XYZ centroid
    if RelType(j)==1
        N_possible_roots = 6;
    else
        N_possible_roots = 2;
    end
        
    %determine the triplet types that apply to the given relationship type
    goodRows = find(TripletTable.Relationship_Type==RelType(j));
    
    %init 
    pR1R2{j} = zeros(N_possible_roots,6);
    for i=1:N_possible_roots
        R1_temp = X(j)*TripletTable.XY_R1(goodRows(i),1) + Y(j)*TripletTable.XY_R1(goodRows(i),2);
        R2_temp = X(j)*TripletTable.XY_R2(goodRows(i),1) + Y(j)*TripletTable.XY_R2(goodRows(i),2);
        pR1R2{j}(i,:) = [R1_temp,R2_temp,j,RelType(j),goodRows(i),Centroid_mag(j)];
    end
end

% Define the distance between clusters to be the Euclidian distance for
% the possible roots that give the lowest value.
D = NaN(N_clusters); %distance
trp_type = zeros(N_clusters); %triplet type for min distance

for i=1:N_clusters    
    for j=1:N_clusters       
            pR1R2_i = pR1R2{i};
            pR1R2_j = pR1R2{j};
            % check all root combinations for the shortest distance
            d_best_temp = Inf;
            trp_type_temp = 1;
            for ii = 1:size(pR1R2_i,1)
                for jj = 1:size(pR1R2_j,1)
                   
                   d_temp = sqrt((pR1R2_i(ii,1)-pR1R2_j(jj,1)).^2 + (pR1R2_i(ii,2)-pR1R2_j(jj,2)).^2);
                   if d_temp < d_best_temp
                       d_best_temp = d_temp;
                       trp_type_temp = ii;
                   end
                   
                end
            end
            D(i,j) = d_best_temp;
            trp_type(i,j) = trp_type_temp;                    
    end
end

D2 = D;
D2(D<max_dist)=1;
D2(D>=max_dist)=0;
%good_trp_types = trp_type(:,best_row);
good_trp_types = mode(trp_type,2);

% now, fix cluster types and compute the distance again:
D3 = Inf(N_clusters); %distance
for i=1:N_clusters    
    for j=1:N_clusters       
            pR1R2_i = pR1R2{i}(good_trp_types(i),:);
            pR1R2_j = pR1R2{j}(good_trp_types(j),:);
            D3(i,j) = sqrt((pR1R2_i(1,1)-pR1R2_j(1,1)).^2 + (pR1R2_i(1,2)-pR1R2_j(1,2)).^2);       
    end
end
D4 = D3;
D4(D3<max_dist)=1;
D4(D3>=max_dist) = 0;

[~,sorted_rows] = sort(sum(D4),'descend'); 
best_row = sorted_rows(1); 
cluster_mask = logical(D4(best_row,:)); %these are the specific clusters that have maximal agreement on root identity
N_good_clusters = sum(cluster_mask);

R1_ = zeros(N_good_clusters,1);
R2_ = zeros(N_good_clusters,1);
Ramp_ = zeros(N_good_clusters,1);

R1 = [];
R2 = [];
Ramp = [];

if sum(cluster_mask)>=minimum_agreeing_clusters
    for i=1:N_clusters
        R1_(i,1) = pR1R2{i}(good_trp_types(i),1);
        R2_(i,1) = pR1R2{i}(good_trp_types(i),2);
        Ramp_(i,1) = pR1R2{i}(good_trp_types(i),6);
    end
    R1 = R1_(cluster_mask);
    R2 = R2_(cluster_mask);
    Ramp = Ramp_(cluster_mask);
end
% R1 =R1_;
% R2 = R2_;
% Ramp = Ramp_;
%try taking all putative roots with sufficient agreement

end
