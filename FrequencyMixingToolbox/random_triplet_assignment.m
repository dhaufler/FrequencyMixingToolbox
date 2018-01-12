function rand_trips = random_triplet_assignment(pR1R2_in)
%...
N_clusters = length(pR1R2_in);
rand_trips = zeros(N_clusters,1);
for i=1:N_clusters
    rand_trips(i,1) = randsample(size(pR1R2_in{i},1),1);
end
end