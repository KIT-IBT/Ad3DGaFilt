% -------------------------------------------------------
%
%    Testscript for the useage of the Ad3DGaFilt filter for optical data
%
%    Ver. 1.0
%
%    Created:       Gerald Schwaderlapp (20.2.2018)
%    Last modified: Gerald Schwaderlapp (20.2.2018)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology (KIT)
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2018 - All rights reserved.
%
% ------------------------------------------------------
%
% Choose data set to filter and run the script:  

clear all; 
close all; 


%% Example (synthetic) optical data
% Load ideal (noise free) data: data_ideal
load('data_ideal.mat')
    
% Load noisy data with -10dB AWGN: data_noisy
load('data_noisy1.mat')
load('data_noisy2.mat')

data_noisy = zeros(82,82,739); 
data_noisy(:,:,1:369) = data_noisy1;
data_noisy(:,:,370:end) = data_noisy2; 

clear data_noisy1 data_noisy2

% Load noisy data -10dB AWGN and baseline: data_noisy_baseline
load('data_noisy_baseline1.mat')
load('data_noisy_baseline2.mat')

data_noisy_baseline = zeros(82,82,739); 
data_noisy_baseline(:,:,1:369) = data_noisy_baseline1;
data_noisy_baseline(:,:,370:end) = data_noisy_baseline2; 

clear data_noisy_baseline1 data_noisy_baseline2

data = data_noisy_baseline;

%% Configuration

% Samplerate of artifical camera
Fs = 868;
% Resolution of artifical camera
res = 16*10^-6; %mm
% 1/Optical magnification
mag = 4;

% time vector
t = (0:1:size(data,3) - 1) / Fs;

%% Filter Settings

% Cutoff frequency for baseline removal
baseline_cutoff = 2; %Hz

% Window length for Welch's method
welchWindow = 0.06; %s
% Energy limit for temporal filtering
tempLimit = 0.95;
% Energy limit for spatial filtering
spatLimit = 0.90;

% Cut data for analysis of filter settings, final filtering is performed an all
    % data
% discard border pixel
disc_px = 5;
% discard signal at beginning and end
disc_sig = 0.15; %s

%% Perform Filtering

% Baseline Removal
data_remB = removeBaseline(data,Fs,baseline_cutoff);

% Ad3DGaFilt
[data_filtered,sigma_spatial,sigma_temp]=Ad3DGaFilt(data_remB,Fs,mag,res,disc_sig,disc_px,welchWindow,tempLimit,spatLimit);

% Normalize data after Ad3DGaFilt  
min_data = repmat(max(data_filtered,[],3),[1 1 size(data_filtered,3)]);
diff_data = repmat(max(data_filtered,[],3)-min(data_filtered,[],3),[1 1 size(data_filtered,3)]);
data_norm = (data_filtered-min_data)./(diff_data);

%% Plots

%do plots for element:
px_xcoord = 20;
px_ycoord = 20;

fh = figure('DefaultAxesFontSize', 20);
set(fh,'Position',[100 100 1000 1000]);

% Plot of ideal signal 
subplot(3,1,1);
plot(t, squeeze(data_ideal(px_xcoord, px_ycoord, :)),'LineWidth', 4);
ylabel('Normalized magnitude (a.u.)');
xlabel('Time (s)');
xlim([t(1) t(size(data_ideal,3))]);
title(strcat('Ideal data for element [',num2str(px_xcoord),',',num2str(px_ycoord),']'));

% Plot of noisy signal with baseline wander
subplot(3,1,2);
plot(t, squeeze(data(px_xcoord, px_ycoord, :)),'LineWidth', 4);
ylabel('Normalized magnitude (a.u.)');
xlabel('Time (s)');
xlim([t(1) t(size(data,3))]);
title(strcat('Raw data for element [',num2str(px_xcoord),',',num2str(px_ycoord),']'));

% Plot of normalized filtered signal
subplot(3,1,3);
plot(t, squeeze(data_norm(px_xcoord, px_ycoord, :)),'LineWidth', 4);
ylabel('Normalized magnitude (a.u.)');
xlabel('Time (s)');
xlim([t(1) t(size(data_norm, 3))]);
title(strcat('Filtered data for element [',num2str(px_xcoord),',',num2str(px_ycoord),']'));
    