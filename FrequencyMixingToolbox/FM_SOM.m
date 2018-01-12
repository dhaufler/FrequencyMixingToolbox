function [SOM_DAT, X_sorted] = FM_SOM(X_in)
%Wrapper function to apply self-organizing map (SOM) to data for frequency
%mixing. The data is clustered by the frequency amplitudes.

%   Input: 
%       X:  N-by-P for N samples and P frequencies (this is how my data is formatted)
%   
%   Output: 
%       X_sorted:   X with samples ordered by class number
%       classes:    N-by-1 vector of sample class number
%       n_classes:  number of classes, K
%       n_totals:   K-by-1 vector of total samples in class 1 to K.

%Step 1: Train SOM to get proto-clusters

[~, SCORE] = pca(log(abs(X_in)));

%x = [((abs(X_in)))]'; %SOM wants data as features-by-samples?
x = SCORE(:,1:20)'; %SOM wants data as features-by-samples?
n_training_iters = 100;

% i want to analyze clustered data for cluster sizes >= 500 samples.
if size(X_in,1) <= 5000
    net = selforgmap([4 2],n_training_iters);
elseif size(X_in,1) > 5000 & size(X_in,1) <= 8000
    net = selforgmap([4 3],n_training_iters);
elseif size(X_in,1) > 8000 & size(X_in,1) <= 10000
    net = selforgmap([4 4],n_training_iters);
elseif size(X_in,1) > 10000 & size(X_in,1) <= 15000
    net = selforgmap([4 6],n_training_iters);
elseif size(X_in,1) > 15000
    net = selforgmap([4 6],n_training_iters);    
end

%disable GUI interface
%net.trainParam.showWindow = false;
net = train(net,x);

SOM_DAT.classes = vec2ind(net(x));
SOM_DAT.n_classes = length(unique(SOM_DAT.classes));

SOM_DAT.n_totals = zeros(SOM_DAT.n_classes,1);
for i=1:SOM_DAT.n_classes
    SOM_DAT.n_totals(i,1) = sum(SOM_DAT.classes==i);
end
[~,SOM_DAT.classes_order]=sort(SOM_DAT.classes);
X_sorted = X_in(SOM_DAT.classes_order,:); %large, keep apart

%get proto-cluster "prototype" spectra (centroids)
SOM_DAT.X_proto_mean = zeros(SOM_DAT.n_classes,size(X_in,2));
SOM_DAT.X_proto_std = zeros(SOM_DAT.n_classes,size(X_in,2));
for i=1:SOM_DAT.n_classes    
    SOM_DAT.X_proto_mean(i,:) = mean(log(abs(X_in(SOM_DAT.classes==i,:))));
    SOM_DAT.X_proto_std(i,:) = std(log(abs(X_in(SOM_DAT.classes==i,:))));
end

% %Step 2: Cluster proto-clusters
% SOM_DAT.D = pdist(SOM_DAT.X_proto_mean);
% SOM_DAT.Z = linkage(SOM_DAT.D,'single');
% dendrogram(SOM_DAT.Z);

end

