function [S_full] = FM_simulate_wrapper1(X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%X(3) = X(3)/100;

%common parameters
params.freqs = 20:20:300;
params.total_t = 120;           %duration in seconds            
params.sampling_rate = 1250;
%signal 1 parameters
params.freq1 = [X(1)-20, X(1)+20];
params.amp1 = [1, 1];
params.N_components1 = 1;
params.Q1 = 20;             %Irregularity parameter: number of control points per second

%signal 2 parameters
params.freq2 = [X(2)-5, X(2)+5];
params.amp2 = [1, 1];
params.N_components2 = 1;
params.Q2 = 20;

%signal 3 parameters
params.freq3 = [75 85];
params.amp3 = [0 0];
params.N_components3 = 3;
params.Q3 = 20;

%noise parameters
params.A_noise = 0.2; %coefficient of noise process
params.k_noise = 0.1; %power_law exponent

%mixing parameters
params.A1 = 1; %coefficient of linear
params.A2 = 0.1;%coefficient of quadratic 
params.A3 = 0;%-0.2; %for cubic term

spec_out = FM_simulate(params);
spec_out_ds = spec_out(:,1:50:end);

%below is used in FM_simulate_wrapper1
%Now, compute FM 
load('Triplet_Table_standard.mat');
FMobj = FreqMix(spec_out_ds',params.freqs,Triplet_Table);
FMobj.FM_FFT_S_full;
S_full = FMobj.fft_pII_analysis.S_full;
end

