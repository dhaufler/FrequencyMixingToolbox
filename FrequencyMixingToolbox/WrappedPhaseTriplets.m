function Data = WrappedPhaseTriplets(Data_IN,freqs)
%Characterize distributions of frequency triplets.
%   Inputs:
%       X:      n-by-p phases in radians (-pi to pi), where n is the number
%               of observations, and p is the number of frequencies.
%       f:      vector of p regularly spaced frequencies. ***IMPORTANT*** f(1) should equal
%               f(k+1) - f(k)
%
%   Outputs:

%parameters
n_freqs = length(freqs);
n_samples = size(Data_IN,1);
n_tests = 13;%3*5*5;
% 
% %check that dimensions are consistent
% if size(X,2)~=n_freqs
%     error('Phases should be (samples x frequencies)')
% end
% 
%init output
Data = [];
Data.f = freqs;
Data.n_samples = n_samples;
Data.xyz_key = zeros(n_tests,3);
%Data.phaseCR = zeros(n_freqs,n_freqs,n_freqs,n_tests);
%Data.phaseCRp = zeros(n_freqs,n_freqs,n_freqs,n_tests);

%Create grids of frequency indices for each frequency relationship. The dimension is set 
%by the length of the frequency vector. Note that all values must be
%natural numbers not larger than N_freqs.

% Lets hard code this...
Data.Frequency_Indices = cell(13,1);
Data.Test = [];
Data.Test.circR = cell(13,1);
Data.Test.circRp = cell(13,1);
[I J] = meshgrid(-4*n_freqs:4*n_freqs,-4*n_freqs:4*n_freqs);

Data.Frequency_Indices{1} = get_good_inds(I+J,I,J,n_freqs); %group 1 (blue)
Data.Frequency_Indices{2} = get_good_inds(I+2*J,I,J,n_freqs); %group 2 (red)
Data.Frequency_Indices{3} = get_good_inds(2*I+J,I,J,n_freqs); %group 3 (red)
Data.Frequency_Indices{4} = get_good_inds(2*I-J,I,J,n_freqs); %group 4 (red)
Data.Frequency_Indices{5} = get_good_inds(I,2*I-2*J,J,n_freqs); %group 5 (green)
Data.Frequency_Indices{6} = get_good_inds(I,J,2*I-2*J,n_freqs); %group 6 (green)
Data.Frequency_Indices{7} = get_good_inds(2*I+2*J,I,J,n_freqs); %group 7 (green)
Data.Frequency_Indices{8} = get_good_inds(2*I-2*J,I,J,n_freqs); %group 8 (green)

%these are control tests and should be computed for all the indices above.
Data.Frequency_Indices{9} = get_good_inds(2*I,I,J,n_freqs); %group 9 
Data.Frequency_Indices{10} = get_good_inds(2*J,I,J,n_freqs); %group 10 
Data.Frequency_Indices{11} = get_good_inds(I,2*J,J,n_freqs); %group 11
Data.Frequency_Indices{12} = get_good_inds(I,I,J,n_freqs); %group 12
Data.Frequency_Indices{13} = get_good_inds(I,J,J,n_freqs); %group 13


Data.xyz_key = ...
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
    
for xxx = 1:n_tests
    
    %X,Y,Z coefficients
    i=Data.xyz_key(xxx,1);
    j=Data.xyz_key(xxx,2);
    k=Data.xyz_key(xxx,3);
    
    X_data = Data_IN(:,Data.Frequency_Indices{xxx}(:,1));
    Y_data = Data_IN(:,Data.Frequency_Indices{xxx}(:,2));
    Z_data = Data_IN(:,Data.Frequency_Indices{xxx}(:,3));
    
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_temp =  circ_r(w,[],[],1);
    [phaseCRp_temp,~] =  circ_rtestND(w,1);
    
    %Controling for nearby frequency relationships and cross-frequency
    %(2-1) relationships
    
    i=Data.xyz_key(9,1);
    j=Data.xyz_key(9,2);
    k=Data.xyz_key(9,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl1 =  circ_r(w,[],[],1);
    [phaseCRp_ctl1,~] =  circ_rtestND(w,1);
    
    i=Data.xyz_key(10,1);
    j=Data.xyz_key(10,2);
    k=Data.xyz_key(10,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl2 =  circ_r(w,[],[],1);
    [phaseCRp_ctl2,~] =  circ_rtestND(w,1);
    
    i=Data.xyz_key(11,1);
    j=Data.xyz_key(11,2);
    k=Data.xyz_key(11,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl3 =  circ_r(w,[],[],1);
    [phaseCRp_ctl3,~] =  circ_rtestND(w,1);
    
    i=Data.xyz_key(12,1);
    j=Data.xyz_key(12,2);
    k=Data.xyz_key(12,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl4 =  circ_r(w,[],[],1);
    [phaseCRp_ctl4,~] =  circ_rtestND(w,1);
    
    i=Data.xyz_key(13,1);
    j=Data.xyz_key(13,2);
    k=Data.xyz_key(13,3);
    w = mod([i.*X_data+j.*Y_data + k.*Z_data]+pi,2*pi)-pi;
    phaseCR_ctl5 =  circ_r(w,[],[],1);
    [phaseCRp_ctl5,~] =  circ_rtestND(w,1);
    
    %write results to output 
    Data.Test.circR{xxx,1} = phaseCR_temp;
    Data.Test.circRp{xxx,1} = phaseCRp_temp;
    Data.Control.circR{xxx,1} = [phaseCR_ctl1; phaseCR_ctl2; phaseCR_ctl3; phaseCR_ctl4; phaseCR_ctl5]; 
    Data.Control.circRp{xxx,1} = [phaseCRp_ctl1; phaseCRp_ctl2; phaseCRp_ctl3; phaseCRp_ctl4; phaseCRp_ctl5]; 
end

end

function IndsOut = get_good_inds(X,Y,Z,N_freqs)
temp = [X(:) Y(:) Z(:)];
good_inds1 = ismember(temp(:,1),1:N_freqs) & ismember(temp(:,2),1:N_freqs) & ismember(temp(:,3),1:N_freqs);
good_inds2 = [temp(:,1) > temp(:,2)] & [temp(:,1) > temp(:,3)] & [temp(:,2) > temp(:,3)];
good_inds = good_inds1 & good_inds2;
IndsOut = temp(good_inds,:);
end

