% -------------------------------------------------------
%
%    removeBaseline - remove baseline wander from optical data
%
%    Ver. 1.0
%
%    Created:       Gerald Schwaderlapp (20.2.2018)
%    Last modified: 
%
%    Institute of Biomedical Engineering
%    Universitaet Karlsruhe (TH)
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2018 - All rights reserved.
%
% ------------------------------------------------------
%
% Input:
%        data: optical data to be filtered
%        Fs: samplerate
%        cutoff_frequency: cutoff frequency for the gaussian highpass in Hz
% Output:
%        data: filtered data
%
%
% Example Usage:
% data = removeBaseline(data,868,1.5)
% 
% Revision history:
%  
function data = removeBaseline(data,Fs,cutoff_frequency)
    %% Gaussian highpass
    % define cut off
    sig_t = round(Fs / (2 * pi * cutoff_frequency));
    % 1D Convolution Filter
        %create gaussian 
        g = fspecial('gauss', [6 * sig_t+1 1], sig_t);    
        %Apply lowpass filter
            %convolution of gaussian and time signal for each pixel     
            lowF = zeros(size(data, 1), size(data, 2), size(data, 3));
            for i = 1:size(data, 1)
                for j = 1:size(data, 2)
                    lowF(i,j, :) = conv(squeeze(data(i,j, :)), g, 'same'); 
                end
            end
    % Difference of data and lowpass filtered data
        data = data - lowF;     
    % Remove median
    for i = 1:size(data, 1)
        for j = 1:size(data, 2)
            data(i, j, :) = data(i, j, :) - median(data(i, j, :));
        end
    end
end


