function circR_analysis = FM_circR(obj)
%Characterize distributions of frequency triplets.


obj.circR_analysis.xyz_key = zeros(13,3);

%Create grids of frequency indices for each frequency relationship. The dimension is set 
%by the length of the frequency vector. Note that all values must be
%natural numbers not larger than N_freqs.

% Lets hard code this...
obj.circR_analysis.Frequency_Indices = cell(13,1);
obj.circR_analysis.Test = [];
obj.circR_analysis.Test.circR = cell(13,1);
obj.circR_analysis.Test.circRp = cell(13,1);
[I J] = meshgrid(-4*obj.N_freqs:4*obj.N_freqs,-4*obj.N_freqs:4*obj.N_freqs);

obj.circR_analysis.Frequency_Indices{1} = get_good_inds(I+J,I,J,obj.N_freqs); %group 1 (blue)
obj.circR_analysis.Frequency_Indices{2} = get_good_inds(I+2*J,I,J,obj.N_freqs); %group 2 (red)
obj.circR_analysis.Frequency_Indices{3} = get_good_inds(2*I+J,I,J,obj.N_freqs); %group 3 (red)
obj.circR_analysis.Frequency_Indices{4} = get_good_inds(2*I-J,I,J,obj.N_freqs); %group 4 (red)
obj.circR_analysis.Frequency_Indices{5} = get_good_inds(I,2*I-2*J,J,obj.N_freqs); %group 5 (green)
obj.circR_analysis.Frequency_Indices{6} = get_good_inds(I,J,2*I-2*J,obj.N_freqs); %group 6 (green)
obj.circR_analysis.Frequency_Indices{7} = get_good_inds(2*I+2*J,I,J,obj.N_freqs); %group 7 (green)
obj.circR_analysis.Frequency_Indices{8} = get_good_inds(2*I-2*J,I,J,obj.N_freqs); %group 8 (green)

%these are control tests and should be computed for all the indices above.
obj.circR_analysis.Frequency_Indices{9} = get_good_inds(2*I,I,J,obj.N_freqs); %group 9 
obj.circR_analysis.Frequency_Indices{10} = get_good_inds(2*J,I,J,obj.N_freqs); %group 10 
obj.circR_analysis.Frequency_Indices{11} = get_good_inds(I,2*J,J,obj.N_freqs); %group 11
obj.circR_analysis.Frequency_Indices{12} = get_good_inds(I,I,J,obj.N_freqs); %group 12
obj.circR_analysis.Frequency_Indices{13} = get_good_inds(I,J,J,obj.N_freqs); %group 13


obj.circR_analysis.xyz_key = ...
    [   1,-1,-1; ...%1
        1,-1,-2; ...%2
        1,-2,-1; ...%3
        1,-2, 1; ...%4
        2,-1,-2; ...%5
        2,-2,-1; ...%6
        1,-2,-2; ...%7
        1,-2, 2; ...%8
        1,-2, 0; ...%9
        1, 0,-2; ...%10
        0, 1,-2; ...%11
        1, -1,0; ...%12
        0, 1, -1];  %13
    
for xxx = 1:13
    
    %X,Y,Z coefficients
    i=obj.circR_analysis.xyz_key(xxx,1);
    j=obj.circR_analysis.xyz_key(xxx,2);
    k=obj.circR_analysis.xyz_key(xxx,3);
    
    X_data = angle(obj.data(:,obj.circR_analysis.Frequency_Indices{xxx}(:,1)));
    Y_data = angle(obj.data(:,obj.circR_analysis.Frequency_Indices{xxx}(:,2)));
    Z_data = angle(obj.data(:,obj.circR_analysis.Frequency_Indices{xxx}(:,3)));
    
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_temp =  circ_r(w,[],[],1);    
    [phaseCRp_temp,~] =  circ_rtestND(w,1);
    
    %Controling for nearby frequency relationships and cross-frequency
    %(2-1) relationships
    
    i=obj.circR_analysis.xyz_key(9,1);
    j=obj.circR_analysis.xyz_key(9,2);
    k=obj.circR_analysis.xyz_key(9,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl1 =  circ_r(w,[],[],1);
    [phaseCRp_ctl1,~] =  circ_rtestND(w,1);
    
    i=obj.circR_analysis.xyz_key(10,1);
    j=obj.circR_analysis.xyz_key(10,2);
    k=obj.circR_analysis.xyz_key(10,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl2 =  circ_r(w,[],[],1);
    [phaseCRp_ctl2,~] =  circ_rtestND(w,1);
    
    i=obj.circR_analysis.xyz_key(11,1);
    j=obj.circR_analysis.xyz_key(11,2);
    k=obj.circR_analysis.xyz_key(11,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl3 =  circ_r(w,[],[],1);
    [phaseCRp_ctl3,~] =  circ_rtestND(w,1);
    
    i=obj.circR_analysis.xyz_key(12,1);
    j=obj.circR_analysis.xyz_key(12,2);
    k=obj.circR_analysis.xyz_key(12,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl4 =  circ_r(w,[],[],1);
    [phaseCRp_ctl4,~] =  circ_rtestND(w,1);
    
    i=obj.circR_analysis.xyz_key(13,1);
    j=obj.circR_analysis.xyz_key(13,2);
    k=obj.circR_analysis.xyz_key(13,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl5 =  circ_r(w,[],[],1);
    [phaseCRp_ctl5,~] =  circ_rtestND(w,1);
    
    %write results to output 
    obj.circR_analysis.Test.circR{xxx,1} = phaseCR_temp;
    obj.circR_analysis.Test.circRp{xxx,1} = phaseCRp_temp;
    obj.circR_analysis.Control.circR{xxx,1} = [phaseCR_ctl1; phaseCR_ctl2; phaseCR_ctl3; phaseCR_ctl4; phaseCR_ctl5]; 
    obj.circR_analysis.Control.circRp{xxx,1} = [phaseCRp_ctl1; phaseCRp_ctl2; phaseCRp_ctl3; phaseCRp_ctl4; phaseCRp_ctl5]; 
end

end

function IndsOut = get_good_inds(X,Y,Z,N_freqs)
temp = [X(:) Y(:) Z(:)];
good_inds1 = ismember(temp(:,1),1:N_freqs) & ismember(temp(:,2),1:N_freqs) & ismember(temp(:,3),1:N_freqs);
good_inds2 = [temp(:,1) >= temp(:,2)] & [temp(:,1) >= temp(:,3)] & [temp(:,2) >= temp(:,3)];
good_inds = good_inds1 & good_inds2;
IndsOut = temp(good_inds,:);
end
