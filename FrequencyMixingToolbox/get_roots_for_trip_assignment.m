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