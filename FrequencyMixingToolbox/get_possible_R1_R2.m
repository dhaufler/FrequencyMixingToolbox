function [R1,R2] = get_possible_R1_R2(X,Y,RelType,TripletTable)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if RelType==1
    N_possible_roots = 6;
else
    N_possible_roots = 2;
end

%get the row numbers that apply
goodRows = find(TripletTable.Relationship_Type==RelType);

R1 = NaN(N_possible_roots,1);
R2 = NaN(N_possible_roots,1);

for i=1:length(goodRows)
    R1(i) = X*TripletTable.XY_R1(goodRows(i),1)+ Y*TripletTable.XY_R1(goodRows(i),2);
    R2(i) = X*TripletTable.XY_R2(goodRows(i),1)+ Y*TripletTable.XY_R2(goodRows(i),2);
end


end

