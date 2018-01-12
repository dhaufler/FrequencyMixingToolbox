%% Construct Triplet Table

% TripletTable
% TripletTable.Triplet_Number
% TripletTable.Relationship_Order
% TripletTable.Relationship_Type
% 
% TripletTable.compX %coefficients of R1 and R2 in the X component
% TripletTable.compY %coefficients of R1 and R2 in the Y component
% TripletTable.compZ %coefficients of R1 and R2 in the Z component
% 
% %coefficients of X, Y, and Z in the frequency and phase relationships
% TripletTable.XYZ_relationship
% 
% %coefficients of X and Y in determining Root1 and Root2
% TripletTable.XY_R1
% TripletTable.XY_R2

%% 
TripletTable = [];
component_types = [     1,0; ... R1
                        0,1; ... R2
                        2,0; ... H1
                        0,2; ... H2
                        1,1; ... E1
                        1,-1]; % E2
                    
% component_types = [     1,0; ... R1
%                         0,1; ... R2
%                         2,0; ... H1
%                         0,2; ... H2
%                         1,1; ... E1
%                         1,-1;... E2
%                         2,-1;...
%                         -1,2;...
%                         1,2;... 
%                         2,1;...
%                         -2,1;...
%                         1,-2];

                        
% try all triplet permutations
k=1;
for x = 1:size(component_types,1)
    for y = 1:size(component_types,1)
        for z = 1:size(component_types,1)
            if (x~=y)&(x~=z)&(y~=z)
                [isGood, relR1R2_condition] = FM_validate_triplet(component_types(x,:),component_types(y,:),component_types(z,:));
                %[hasPairwise, relR1R2_condition] = FM_check_pairwise(component_types(x,:),component_types(y,:),component_types(z,:));
                if isGood ==1
                    TripletTable(k).relative_R1_R2_frequency = relR1R2_condition;
                  	TripletTable(k).compX = component_types(x,:);
                 	TripletTable(k).compY = component_types(y,:);
                 	TripletTable(k).compZ = component_types(z,:);
                    
                    TripletTable(k).XYZ_relationship = FM_determine_relationship_type(component_types(x,:),component_types(y,:),component_types(z,:));
                    TripletTable(k).Relationship_Order = sum(abs(TripletTable(k).XYZ_relationship));
                    
                    %Below are the X and Y coefficients to recover R1 and R2 from X and Y.
                    a = component_types(x,1); b = component_types(x,2);
                    c = component_types(y,1); d = component_types(y,2);
                    [R1_X, R1_Y,R2_X,R2_Y] = FM_solve_roots(component_types(x,:),component_types(y,:));
                    TripletTable(k).XY_R1 = [R1_X,R1_Y];
                    TripletTable(k).XY_R2 = [R2_X,R2_Y];
                    k=k+1
                end
            end
        end
    end
end

Triplet_Table_ = struct2table(TripletTable);
Triplet_Table_ = sortrows(Triplet_Table_,'Relationship_Order','ascend');

%determine the number of distinct XYZ_relationships
[distinct_XYZ, IA,IC] = unique(Triplet_Table_.XYZ_relationship(:,:),'rows','stable');
%IC assigns a relationship type number to each triplet, but these are
%arbitrary labels. I want labels increasing with relationship order.
n_distinct_XYZ = size(distinct_XYZ,1);
for i=1:size(Triplet_Table_,1)
    Triplet_Table_.Relationship_Type(i) = IC(i);
end
Triplet_Table_ = sortrows(Triplet_Table_,'Relationship_Type','ascend');
for i=1:size(Triplet_Table_,1)
    Triplet_Table_.Triplet_Number(i) = i;
end
Triplet_Table = Triplet_Table_(:,{'Triplet_Number','Relationship_Type','Relationship_Order','relative_R1_R2_frequency',...
    'compX','compY','compZ','XYZ_relationship','XY_R1','XY_R2'})

%% scrap

distinct_XYZ = unique(Triplet_Table_.XYZ_relationship(:,:),'rows','stable');
[R1_X, R1_Y,R2_X,R2_Y] = FM_solve_roots([2 0],[1 1])
 [isGood, relR1R2_condition] = FM_validate_triplet([1 0],[0 1],[1,-1]);
