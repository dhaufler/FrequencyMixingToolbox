function [X_out,leafOrder] = FM_group_by_power(X)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

X_abs = (abs(X));
D = pdist(X_abs);
tree = linkage(D,'average');
leafOrder = optimalleaforder(tree,D);
X_out = X(leafOrder,:);
end

