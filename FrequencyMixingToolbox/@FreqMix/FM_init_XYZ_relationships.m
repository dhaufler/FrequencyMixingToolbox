function obj = FM_init_XYZ_relationships(obj)
%UNTITLED Summary of this function goes here

%Create grids of frequency indices for each frequency relationship. The dimension is set 
%by the length of the frequency vector. Note that all values must be
%natural numbers not larger than N_freqs.

obj.Frequency_Indices = cell(obj.N_distinct_XYZ_relationships,1);
[I J] = meshgrid(-6*obj.N_freqs:6*obj.N_freqs,-6*obj.N_freqs:6*obj.N_freqs);
for i=1:obj.N_distinct_XYZ_relationships
    a = obj.distinct_XYZ_relationships(i,1);
    b = obj.distinct_XYZ_relationships(i,2);
    c = obj.distinct_XYZ_relationships(i,3);
    obj.Frequency_Indices{i} = get_good_inds(I,J,(I*(-a)+J*(-b))/c,obj.N_freqs); %group 1 (blue)
end


end

function IndsOut = get_good_inds(X,Y,Z,N_freqs)
temp = [X(:) Y(:) Z(:)];
good_inds1 = ismember(temp(:,1),1:N_freqs) & ismember(temp(:,2),1:N_freqs) & ismember(temp(:,3),1:N_freqs);
good_inds2 = [temp(:,1) >= temp(:,2)] & [temp(:,1) >= temp(:,3)] & [temp(:,2) >= temp(:,3)];
good_inds = good_inds1 & good_inds2;
IndsOut = temp(good_inds,:);
end