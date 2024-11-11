function [comp0,COMP] = DecomposeFourier_func(x, fs, f1)

overlapRatio   = 0.95; %  Put it back to 0.99
windowFunction = @hann; % Window function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('sampledata10.mat'); % Contains: f1, fs, x(t) [Pa] data
t  = [0:length(x)-1]/fs;
oscillatingSignal = x';

% Sliding window analysis
T1 = 1 / f1;        % Period of the fundamental frequency for component 1
windowMultiple = 2; % Def (2 or 4)
windowLength   = ceil(fs * T1 * windowMultiple);

hopSize = ceil((1 - overlapRatio) * windowLength); % Hop size  (round-->ceil, 2024-03-04)
numWindows = floor((length(oscillatingSignal) - windowLength) / hopSize) + 1;

% Initialize variables to store the results
amplitudeData1 = zeros(numWindows, 1);
amplitudeData2 = zeros(numWindows, 1);
amplitudeData3 = zeros(numWindows, 1);
timeVector = zeros(1, numWindows);
MXFFT = zeros(windowLength,numWindows);


% Perform the sliding window analysis
for i = 1:numWindows
    startIndex = (i - 1) * hopSize + 1;
    endIndex = startIndex + windowLength - 1;
    % Ensure endIndex is within bounds
    endIndex = min(endIndex, length(oscillatingSignal));
    
    windowedSegment = oscillatingSignal(startIndex:endIndex) .* windowFunction(windowLength);
      
    MXFFT(:,i) = fft(windowedSegment)/sqrt(windowLength);
      
    % Calculate the corresponding time value as the midpoint of the window
    timeVector(i) = mean(t(startIndex:endIndex));
end


% RECOMPOSE SIGNAL IN iFFT %++++++++++++++++++++++++++++++++++++++++
Nbins  = min(size(MXFFT));
Ncomps = 5;
COMP   = zeros(length(oscillatingSignal),Ncomps);

for idx_comp = 1 : Ncomps

mask                                     = ones(Nbins,1);
mask(round(idx_comp*f1*windowLength/fs))         = 0;
mask(end - round(idx_comp*f1*windowLength/fs) +1) = 0;
mask                                = repmat(mask, [1,max(size(MXFFT))]);
MXFFTmask             = MXFFT;
MXFFTmask(find(mask)) = [0];
comp = zeros(length(oscillatingSignal),1);

    for i = 1:numWindows
        startIndex = (i - 1) * hopSize + 1;
        endIndex = startIndex + windowLength - 1;
        endIndex = min(endIndex, length(comp));

        tmp = real(ifft(MXFFTmask(:,i)))*sqrt(windowLength);
        tmp = tmp(:);
        comp(startIndex:endIndex) = comp(startIndex:endIndex) + tmp;
    end     
    COMP(:,idx_comp) = comp;
end    

% Comp 0 (DC)
mask= ones(Nbins,1);
mask(1)   = 0;
mask(end) = 0;
mask = repmat(mask, [1,max(size(MXFFT))]);
MXFFTmask = MXFFT;
MXFFTmask(find(mask)) = [0];
comp0 = zeros(length(oscillatingSignal),1);

for i = 1:numWindows
    startIndex = (i - 1) * hopSize + 1;
    endIndex = startIndex + windowLength - 1;    
    endIndex = min(endIndex, length(comp0));
        
    tmp = real(ifft(MXFFTmask(:,i)))*sqrt(windowLength);
    tmp = tmp(:);    
    comp0(startIndex:endIndex) = comp0(startIndex:endIndex) + tmp;    
end

% SCALE THEM ============
TOTAL = comp0 + sum(COMP,2);
fac = max(abs(x))/max(abs(TOTAL));
comp0 = comp0*fac;
COMP  = COMP*fac;
TOTAL = TOTAL*fac;
