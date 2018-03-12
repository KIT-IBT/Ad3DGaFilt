% -------------------------------------------------------
%
%    Ad3DGaFilt - adaptive 3D Gauss filter for optical data
%
%    Ver. 1.0
%
%    Created:           Gerald Schwaderlapp (20.2.2018)
%    Last modified:     Gerald Schwaderlapp (20.2.2018)
%
%    Institute of Biomedical Engineering (IBT)
%    Karlsruhe Institute of Technology (KIT)
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2018 - All rights reserved.
%
% ------------------------------------------------------
%
% [data,sigma_spatial,sigma_temporal]= Ad3DGaFilt(data,Fs,mag,res,disc_sig,disc_px,welchWindow,tempLimit,spatLimit)
%
% Input:
%        data: optical data to be filtered
%        Fs: samplerate
%        mag: optical magnification
%        res: resolution in m
%        disc_sig: Part of the signal discarded at the beginning and the
%                  end for the analysis
%        disc_px: Number of pixels discarded at the border for analysis
%        welchWindow: Window size for Welch's method
%        tempLimit: Energy limit for temporal filtering
%        spatLimit: Energy limit for spatial filtering
% Output:
%        data: filtered data
%        sigma_spatial: spatial sigma used for the filter kernel
%        sigma_temporal: temporal sigma used for the filter kernel
%
% 
% Example Usage:
% [data,sigma_spatial,sigma_temporal]=Ad3DGaFilt(data,868,4,16*10^-6,0.15,5,0.06,0.95,0.95)
% 
% Revision history:
%  

function [data,sigma_spatial,sigma_temporal]=Ad3DGaFilt(data,Fs,mag,res,disc_sig,disc_px,welchWindow,tempLimit,spatLimit)

%% Discard data for analysis
    starttime = round(Fs * disc_sig);
    if size(data, 3) - starttime >= 2 * starttime
        data_disc = data(disc_px:end - disc_px,disc_px:end - disc_px,starttime:end - starttime);
    else
        error('Failed to cut data for analysis!')
    end
    
%% Temporal analysis
    %calculate spatial mean value
    spatial_mean = zeros(size(data_disc, 3), 1);
    for i=1:size(data_disc, 3)
        temp = data_disc(:, :, i);
        spatial_mean(i) = mean(temp(:));
    end
    %time window for Welch's method
    window = round(Fs * welchWindow);
    %overlap
    overlap = round(0.5 * window);
    %zero-padding
    nfft = 8 * 2^nextpow2(window);
    %Welch's method
    [welchPxx,ftemporal] = pwelch(spatial_mean, window, overlap, nfft, Fs, 'power');
    %calculate frequency for tempLimit% of energy
    i=1;
    while sum(welchPxx(1:i))/sum(welchPxx(:)) < tempLimit
        temporal = ftemporal(i);
        i = i + 1;
    end
    %calculate temporal sigma
    sigma_temporal = round(Fs / (temporal * 2 * pi));
    %point in time for maximal wave spread in the image
    maxT = find(spatial_mean == min(spatial_mean));

%% Spatial analysis
    % temporal gauss filtering
        % create gaussian kernel for temporal filtering
        tgauss = fspecial('gauss',[6*round(sigma_temporal)+1 1],round(sigma_temporal));
    % convolution of gaussian and time signal for each element
    for i=1:size(data_disc,1)
        for j=1:size(data_disc,2)
            data_disc(i,j,:) = conv(squeeze(data_disc(i, j, :)), tgauss,'same');
        end
    end   
    % ACF of image at maxT
    [data_disc_x,data_disc_y] = size(data_disc(:, :, maxT(1)));
    cor = xcorr2(data_disc(:, :, maxT(1)));  
    data_disc_xcorr = cor(ceil(data_disc_x / 2):floor(3 * data_disc_x / 2) - 1,...
        ceil(data_disc_y / 2):floor(3 * data_disc_y/2) - 1);
    % properties for fft and transformation in polar coordinates
    N = 8 * 2^nextpow2(size(data_disc_xcorr, 1));
    hN = N / 2;
    [X, Y] = meshgrid(-hN:hN - 1,-hN:hN - 1);
    [~, rho] = cart2pol(X, Y);
    rho = round(rho);
    radAvg = zeros(hN, 1);
    % spatial radial frequency vector 
    fspatial = (1 / (res*mag)) * (0:hN - 1)/N;
    % 2D-FFT of autocorrelation
    mFFT = fft2(data_disc_xcorr, N, N);
    mFFTabs = fftshift((abs(mFFT / (N * N))));
    % radial averaged power spectrum
    for r=0:hN-1
        radAvg(r+1, 1) = nanmean(mFFTabs(rho == r));
    end
    % calculate frequency for spatLimit% of energy 
    i = 1;
    while sum(radAvg(1:i)) / sum(radAvg(:)) <= spatLimit
        spatial = fspatial(i);
        i = i + 1;
    end
    
    %calculate spatial sigma
    sigma_spatial=floor((1 / (res * mag)) / (spatial * 2 * pi));

%% Apply final filtering
    % kernel
    sigma = [abs(sigma_spatial), abs(sigma_spatial), abs(sigma_temporal)]; 
    % filter
    data = imgaussfilt3(data,sigma,'Filtersize',6*ceil(sigma) + 1);

end    