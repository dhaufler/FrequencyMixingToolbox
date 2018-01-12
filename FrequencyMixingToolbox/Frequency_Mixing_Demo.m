%% FREQUENCY MIXING EXAMPLE SCRIPT. SEPTEMBER 12, 2017. DARRELL HAUFLER

%% Simulate mixing process of two signals and get spectrogram

%common parameters
freqs = 2:2:320;
total_t = 30;           %duration in seconds            
sampling_rate = 1250;

%signal 1 parameters
freq1 = [139 141];
amp1 = [0.8 1.2];
N_components1 = 1;
Q1 = 5;             %Irregularity parameter: number of control points per second

%signal 2 parameters
freq2 = [90 110];   %Fig. 2e: 85-115 Fig.2b 30-50
amp2 = [0.0 0.0];     %Fig. 2e: 0 0.5 Fig.2b 0-0.25
N_components2 = 1;  %Fig. 2e: 1 Fig.2b 1
Q2 = 1;            %Fig. 2e: 20 Fig.2b 20

%signal 3 parameters
freq3 = [39 41];
amp3 = [0.8 1.2];
N_components3 = 1;
Q3 = 5;

%noise parameters
A_noise = 0.2; %coefficient of noise process
k_noise = 0.1; %power_law exponent

%mixing parameters
A1 = 1; %coefficient of linear      %Fig.2b 1
A2 = 5;%coefficient of quadratic %Fig.2b 0.25
A3 = 0;%-0.2; %for cubic term

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
%additive_noise =0.2*randn(length(root1),1);
noiseObj = dsp.ColoredNoise('SamplesPerFrame',length(params.t), 'InverseFrequencyPower',k_noise);
additive_noise = A_noise*noiseObj(); %instantiation of noise

%linear_sum      = A1*(root1' + root2')+additive_noise;
nonlinear_sum   = A1*(root1' + root2'+root3')+A2*(root1' + root2'+root3').^2 +A3*(root1' + root2' +root3').^3+additive_noise;
%linear_sum       = A1*(root1' + root2'+root3')+A2*(root1').^2 +A2*(root2').^2  +additive_noise;

%get spectrogram
%X_l     = WaveletSpec(linear_sum,1/1250,freqs,'WAVELETSPEC',{'morlet' 15},'Normalize');
X_nl    = WaveletSpec(nonlinear_sum,1/1250,freqs,'WAVELETSPEC',{'morlet' 50},'Normalize');

%%plot spectrograms of linear vs. nonlinearly combined signals
%figure; 
%S_max = max(abs([X_l.Spectrum(:);X_nl.Spectrum(:)]));
%subplot(1,2,1); imagesc(params.t,freqs,abs(X_nl.Spectrum),[0 S_max]); axis xy
%subplot(1,2,2); imagesc(params.t,freqs,abs(X_l.Spectrum),[0 S_max]); axis xy

%plot just nl
figure; 
imagesc(params.t,freqs,abs(X_nl.Spectrum)); axis xy
ylim([0 240])
title('Cubic')
set(gca,'XTick',[0 0.5 1 1.5 2]);
set(gca,'YTick', [0 50 100 150 200])

%%  Analyize simulated data

% Step 1: Create Frequency Mixing (FM) object
X = X_l.Spectrum(:,1:200:end)'; %spectral data, samples-by-freqs
%X = channelSamples(:,:)';
%tic;myFMquad_good3 = FreqMix(X(:,1:2:end),1:2:480,Quad_Table);toc
tic;myFM4 = FreqMix(X(:,4:4:end),4:4:320,Triplet_Table);toc
%tic;myFMcubic5 = FreqMix(X(:,1:2:end),1:2:480,Cubic_Table);toc
myFMnc = FreqMix(X,1:1:320);
myFMs = FreqMix(X_s,1:1:320);

myFMquad3.image_R1_R2_by_triplets(1:20);
myFMquad_good3.image_R1_R2_by_triplets(1:20);
myFMcubic5.image_R1_R2_by_triplets(45);

type1 = [1 3 5 6 7 10 11 12 15 16 18 20];
Only_type1 = [7 10 11 12 15 16 18 20];
only_unique = [];
%%  Plot 3D scatter of phase distribution (this is used in FIGURE_S1_Nov22_work
figure;
sample_inds = randsample(size(X_nl.Spectrum,2),1000);
phases_X = angle(X_nl.Spectrum(20,sample_inds))'; %40 Hz 
phases_Y = angle(X_nl.Spectrum(50,sample_inds))'; %100 Hz 
phases_Z = angle(X_nl.Spectrum(70,sample_inds))'; %140 Hz 
scatter3(phases_X,phases_Y,phases_Z,20,'k','filled','MarkerFaceAlpha',0.35,'MarkerEdgeAlpha',0); 
xlim([-pi pi])
ylim([-pi pi])
zlim([-pi pi])
xlabel('X-phase')
ylabel('Y-phase')
zlabel('Z-phase')
set(gca,'XTick',[-pi -pi/2 0 pi/2 pi])
set(gca,'YTick',[-pi -pi/2 0 pi/2 pi])
set(gca,'ZTick',[-pi -pi/2 0 pi/2 pi])

daspect([1 1 1])
box on; grid on
view(60,30)
%% Compare the above result to surrogate data that preserves the pairwise relationships
X = X_(inds_order==i,:);
X_s = FM_surrogate([X]);
myFMs = FreqMix(X_s,f);
myFMs.image_R1_R2_by_triplets(1:20)

%% loaded 'LFP'
t = 1:1250*2000;
f = 1:1:320;
S = WaveletSpec(LFP(t,3),1/1250,f,'WAVELETSPEC',{'morlet' 15});
X_=S.Spectrum(5:5:end,1:500:end)';
%[X_sorted,inds_order] = FM_SOM((X_(:,1:1:240)));
X = channelSamples';
[SOM_DAT,X_sorted] = FM_SOM(X(:,:));
figure; imagesc(size(X,1),1:320,((((abs(X_sorted)))))'); axis xy;
%
%all_avs = zeros(681,681,1);
%rel1_avs = zeros(68,68,1);
for i=1:SOM_DAT.n_classes
%myFM = FreqMix(X_(i*200:(i*200+200),:),f);
%if sum(inds_order==i)>150&sum(inds_order==i)<400
if sum(SOM_DAT.classes==i)>1500
i
dat_temp = X(SOM_DAT.classes==i,:);
dat_temp2 = dat_temp(randsample(size(dat_temp,1),1500),:);
myFM = FreqMix(dat_temp2,1:320,Quad_Table);
%myFMc = FreqMix(X(SOM_DAT.classes==i,:),f,Cubic_Table);
%size(X(sample_class==i,:))
%myFMc.image_R1_R2_by_triplets([53]);
myFM.image_R1_R2_by_triplets([1:20]);
drawnow
end
end
% all_avs(:,:,i) = myFM.circR_analysis.standard_Rvals_mean;
% rel1_avs(:,:,i) = myFM.circR_analysis.IJ_data_no_control(:,:,1);
% figure(199); imagesc(mean(rel1_avs,3))
% figure(155); imagesc(1:f(end)/2,1:f(end/2),mean(all_avs,3),[0 0.2]); axis xy
% drawnow
% waitforbuttonpress
% end

%%
myFM.image_R1_R2_by_triplets([1 3 5 6 7 10 11 12 15 16 18 20])
myFM.image_R1_R2_by_triplets([1 3 4 6 7 8 11 12 13 16 18 20])
%myFMnc.image_R1_R2_by_triplets([2 3 4 6 8 9 11 12 14 17 19 20])

%myFMnc.image_R1_R2_by_triplets([16])

