function [R1,R2,Ramp,R_inds,median_R1_,median_R2_] = get_possible_R1_R2_from_set2(X,Y,RelType,Centroid_mag, TripletTable,max_distance,min_clu_N,N_iterations)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

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

optimal_CV = Inf;
best_trips = [];
for w=1:10
    rand_trip = random_triplet_assignment(pR1R2);
    new_trips = rand_trip;
    local_CV_best = Inf;
    for i=1:N_iterations
        temp_trips = new_trips;
        k = randsample(N_clusters,1);
        n_possible_trips = size(pR1R2{k},1);
        if n_possible_trips==2
            new_k_trip = mod(new_trips(k,1),2)+1;
        else
            new_k_trip = mod(new_trips(k,1),6)+1;
        end     
        temp_trips(k,1) = new_k_trip;
        [R1_temp,R2_temp,~] = get_roots_for_trip_assignment(pR1R2,temp_trips);         
        mean_R1_temp = mean(R1_temp);
        mean_R2_temp = mean(R2_temp);
        length_R1R2 = sqrt(mean_R1_temp.^2+mean_R2_temp.^2);
        R1_temp_subtract_mean = R1_temp - mean_R1_temp;
        R2_temp_subtract_mean = R2_temp - mean_R2_temp;
        %local_CV_temp = sum(R1_temp_subtract_median(R1_temp<160&R2_temp<160).^2 + R2_temp_subtract_median(R1_temp<160&R2_temp<160).^2)/(0.001+sum(R1_temp<160&R2_temp<160)*length_R1R2);
        %local_CV_temp = mean(sqrt(R1_temp_subtract_mean.^2 + R2_temp_subtract_mean.^2))/(length_R1R2);        
        local_CV_temp = mean((R1_temp_subtract_mean.^2 + R2_temp_subtract_mean.^2));
        %local_CV_temp = sum((R1_temp_subtract_mean.^2 + R2_temp_subtract_mean.^2));
        %local_CV_temp = std(R1_temp)/mean(R1_temp) + std(R2_temp)/mean(R2_temp);
        %local_CV_temp = std(R1_temp) + std(R2_temp);
        %total_CV_temp = sum(R1_temp_subtract_mean.^2 + R2_temp_subtract_mean.^2)/(N_clusters*length_R1R2);
        if  local_CV_temp<local_CV_best
            new_trips = temp_trips;
            local_CV_best = local_CV_temp;
        end    
    end
    if local_CV_best<optimal_CV
        optimal_CV = local_CV_best;
        best_trips = new_trips;
    end
end

%get R1 and R2
[R1_,R2_,Ramp_] = get_roots_for_trip_assignment(pR1R2,best_trips);
% only keep within max_distance from the median R1 and R2
median_R1_ = median(R1_);
median_R2_ = median(R2_);
R_dists = zeros(N_clusters,1);
for i=1:N_clusters
    R_dists(i) = sqrt((R1_(i)-median_R1_).^2+(R2_(i)-median_R2_).^2);
end
goodInds = R_dists< max_distance;

if sum(goodInds)>=min_clu_N
    R1 = R1_(goodInds,1);
    R2 = R2_(goodInds,1);
    Ramp = Ramp_(goodInds,1);
    R_inds = goodInds;
else
    R1=[];
    R2=[];
    Ramp=[];
    R_inds = logical(zeros(N_clusters,1));
end
end


function rand_trips = random_triplet_assignment(pR1R2_in)
%...
N_clusters = length(pR1R2_in);
rand_trips = zeros(N_clusters,1);
for i=1:N_clusters
    rand_trips(i,1) = randsample(size(pR1R2_in{i},1),1);
end
end

function [R1,R2,Ramp] = get_roots_for_trip_assignment(pR1R2_in,trps_in)
%...
N_clusters = length(pR1R2_in);    
R1 = zeros(N_clusters,1);
R2 = zeros(N_clusters,1);
Ramp = zeros(N_clusters,1);

for i=1:N_clusters
    R1(i,1) = pR1R2_in{i}(trps_in(i),1);
    R2(i,1) = pR1R2_in{i}(trps_in(i),2);
    Ramp(i,1) = pR1R2_in{i}(trps_in(i),6);
end
end
