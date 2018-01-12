function [spec_out] = FM_simulate(P)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%common parameters
freqs = P.freqs;
total_t = P.total_t;           %duration in seconds            
sampling_rate = P.sampling_rate;

%signal 1 parameters
freq1 = P.freq1;
amp1 = P.amp1;
N_components1 = P.N_components1;
Q1 = P.Q1;             %Irregularity parameter: number of control points per second

%signal 2 parameters
freq2 = P.freq2;
amp2 = P.amp2;
N_components2 = P.N_components2;
Q2 = P.Q2;

%signal 4 parameters
freq3 = P.freq3;
amp3 = P.amp3;
N_components3 = P.N_components3;
Q3 = P.Q3;

%noise parameters
A_noise = P.A_noise; %coefficient of noise process
k_noise = P.k_noise; %power_law exponent

%mixing parameters
A1 = P.A1; %coefficient of linear
A2 = P.A2 ;%coefficient of quadratic 
A3 = P.A3;%-0.2; %for cubic term

%###############
%Construct signal 1
params.Fs = sampling_rate;
params.dt = 1/params.Fs;
params.t = params.dt:params.dt:total_t;
params.ctl_pts_t    = linspace(params.t(1),params.t(end),round(Q1*total_t));
params.IF_low       = ones(size(params.ctl_pts_t))*freq1(1);
params.IF_high      = ones(size(params.ctl_pts_t))*freq1(2);
params.IA_low       = ones(size(params.ctl_pts_t))*amp1(1);
params.IA_high      = ones(size(params.ctl_pts_t))*amp1(2);
params.n_signals = N_components1;
root1_ = my_random_signal(params);
root1 = sum(root1_,1);

%construct signal 2
params2.Fs = sampling_rate;
params2.dt = 1/params2.Fs;
params2.t = params2.dt:params2.dt:total_t;
params2.ctl_pts_t    = linspace(params2.t(1),params2.t(end),round(Q2*total_t));
params2.IF_low       = ones(size(params2.ctl_pts_t))*freq2(1);
params2.IF_high      = ones(size(params2.ctl_pts_t))*freq2(2);
params2.IA_low       = ones(size(params2.ctl_pts_t))*amp2(1);
params2.IA_high      = ones(size(params2.ctl_pts_t))*amp2(2);
params2.n_signals = N_components2;
root2_ = my_random_signal(params2); 
root2 = sum(root2_,1);

%construct signal 3
params3.Fs = sampling_rate;
params3.dt = 1/params3.Fs;
params3.t = params3.dt:params3.dt:total_t;
params3.ctl_pts_t    = linspace(params3.t(1),params3.t(end),round(Q3*total_t));
params3.IF_low       = ones(size(params3.ctl_pts_t))*freq3(1);
params3.IF_high      = ones(size(params3.ctl_pts_t))*freq3(2);
params3.IA_low       = ones(size(params3.ctl_pts_t))*amp3(1);
params3.IA_high      = ones(size(params3.ctl_pts_t))*amp3(2);
params3.n_signals = N_components3;
root3_ = my_random_signal(params3); 
root3 = sum(root3_,1);

%add noise using built-in colored noise generator
noiseObj = dsp.ColoredNoise('SamplesPerFrame',length(params.t), 'InverseFrequencyPower',k_noise);
additive_noise = A_noise*noiseObj(); %instantiation of noise

nonlinear_sum   = A1*(root1' + root2'+root3')+A2*(root1' + root2'+root3').^2 +A3*(root1' + root2' +root3').^3+additive_noise;

%get spectrogram
X_nl    = WaveletSpec(nonlinear_sum,1/1250,freqs,'WAVELETSPEC',{'morlet' 15},'Normalize');
spec_out = X_nl.Spectrum;
end

