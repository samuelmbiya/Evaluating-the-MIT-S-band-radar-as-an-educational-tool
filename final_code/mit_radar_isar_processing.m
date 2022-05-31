function [] = mit_radar_isar_processing(wavFile)

% Author: Samuel Mbiya
%
% Contributions: Gregory Charvat wrote the original code for signal processing and the generation
% of range vs time plots for the Lincoln Laboratory MIT IAP Radar Course
%
% function [] = mit_radar_isar(wavFile)
%
% Produces multiple ISAR images using radar data collected from the MIT Coffeee can radar
%
% Parameters for ISAR imaging demonstration:
%
% turn_table_speed = rotating speed of the rotating platform in RPM
% no_of_pulses_to_consider = number of pulses in the range profiles to consider (N = Dwell_time/PRI)
% fStart = LFM starting frequency (Hz)
% fStop = LFM stopping frequency (H)
% maxRange = maximum range to display in ISAR image
no_pulses_to_consider = 16;
turn_table_speed = 3.635;

fStart = 2350e6;%Hz) LFM start frequency 
fStop = 2530e6;%(Hz) LFM stop frequency
maxRange = 5;%(m) maximum range to display in the ISAR image

% constants

c = 299e6; % (m/s) speed of light
Tp = 10e-3; % (s) minimum pulse length

numPad = 64; % number of samples to pad for bandlimited interpolation & shifting
             
ovsTrg = 16; % oversampling factor applied when interpolating the trigger signal

ovsRng = 1; % oversampling factor applied when interpolating the range data 

nPulseCancel = 3; % number of pulses to use for canceller 
     

% Read the raw wave data
fprintf('Using %g MHz bandwidth\n', (fStop-fStart)*1e-6);
fprintf('Loading WAV file...\n');
[Y,Fs] = audioread(wavFile,'native');

%% Derived parameters
Np = round(Tp * Fs); % Number of samples in the chirp signal (or per pulse)
BW = fStop - fStart; % (Hz) transmit bandwidth
delta_r = c/(2*BW);  % (m) range resolution

trig = -Y(:,1); % the trigger signal on the first channel 
s = -Y(:,2); % the raw mixer output on the second channel
clear Y;

% Estimate the actual chirp pulse length from the measured data and store in variable Tp 

% Estimate the index of the rising edge of the trigger and store in variable pulseStarts
% parse the trigger signal (look for threshold crossings)

fprintf('Parsing the recording...\n');
pulseTrig = (trig > 0);

pulseSum = sum(pulseTrig(1:Np));

pulseStarts = []; %pulseStarts contains the index of the rising edge in samples

% Find the positions of the rising edges in sync pulse
for ii = Np+1:numel(trig)-Np-1        
    if (pulseTrig(ii) && pulseSum==0)
        pulseStarts = [pulseStarts; ii];
    end
   
    % update the running sum
    pulseSum = pulseSum + pulseTrig(ii) - pulseTrig(ii - Np);
end
clear pulseTrig;

% Set the pulse width
Np = round(min(diff(pulseStarts))/2); % Np = number of samples of the chirp pulse

Tp = Np / Fs;
fprintf('Measured pulse width of %g ms \n', Tp*1e3);

%% Pre-compute some windows and other vectors 
Nrange = floor(ovsRng*Np/2); % number of output range samples
dataRange = (0:Nrange-1).' * (delta_r/ovsRng); % labelling the range axes in meters
dataRange = dataRange(dataRange <= maxRange); % labelling the range axes in meters upto maxRange in meters
Nrange_keep = numel(dataRange); % number of range bins which are less than maxRange in meters
%
% Setup window functions 
rngWin = hann_window(Np); % the window applied to reduce range sidelobes
padWin = sin( (1:numPad).'/(numPad+1) * pi/2) .^2; % the window applied to the padded data
trgWin = hann_window(numPad*2+1); % the window applied to the trigger data

%  Obtain the number of pulses in the data
nSamples = numel(s);
pulseStarts = pulseStarts(pulseStarts+Np+numPad <= nSamples);
numPulses = numel(pulseStarts); 
fprintf('Found %d pulses\n',numPulses);

% process pulses into a data matrix
sif = zeros(Nrange_keep,numPulses); % sif - Range lines in matrix form 

fprintf('Processing pulse data...\n');

% Generate range profiles 
for pIdx = 1:numPulses
    % bandlimited interpolate the trigger signal
	tmp = double(trig(pulseStarts(pIdx) + (-numPad:numPad))) .* trgWin; 
                 
    interpTmp = fft_interp(tmp,ovsTrg);
    interpTmp = interpTmp( (numPad*ovsTrg + 1) + (-ovsTrg:ovsTrg) );
    interpOffs = (-ovsTrg:ovsTrg)/ovsTrg;
    myIdx = find(diff(sign(interpTmp))==2)+1;
    tmp2 = interpTmp( myIdx + (-1:0) );
	
    % linear interpolate to find the zero crossing
    fracOffset = -(interpOffs(myIdx) - tmp2(2)/(tmp2(2)-tmp2(1)) / ovsTrg);
    
    % time-align the data to the trigger event (the zero crossing) 
    % Apply non-integer time-shift in the frequency domain by multiplying by a phase ramp signal
    cInds = pulseStarts(pIdx) + (-numPad:(Np+numPad-1)); % pulseStarts(pIdx) = contains estimate of indx of zero crossing
    tmp = double(s(cInds)); %tmp = data received after chirp pulse
	
    %get samples from data from (rising edge-NumPad: NumSamplesChirpPulse + NumPad) 
    tmp(1:numPad) = tmp(1:numPad) .* padWin;
    tmp(end:-1:(end-numPad+1)) = tmp(end:-1:(end-numPad+1)) .* padWin;
   
    % time delay applied in the frequency domain below
    tmp = fft(tmp); 
    tmp = tmp .* exp( -1j*(0:(Np+2*numPad-1)).'/(Np+2*numPad)*2*pi*fracOffset );
    tmp = ifft(tmp,'symmetric');
    
    % compute & scale range data from the time-aligned mixer output
    tmp = ifft(tmp(numPad + (1:Np)) .* rngWin, 2*Nrange); %Perform IFFT to get range line
    sif(:,pIdx) = tmp(1:Nrange_keep); % sif = range profiles in matrix form
end

clear s trig;
sif = sif.'; 

% Display the range profiles

% display the RTI
figure;
imagesc(dataRange,(0:numPulses-1)*Tp*2,20*log10(abs(sif)));
ylabel('Time (s)');
xlabel('Range (m)');
title('RTI without clutter rejection');
colormap(jet(256));
colorbar;
axis xy;

% apply the N-pulse canceller
mti_filter = -ones(nPulseCancel,1)/nPulseCancel;
midIdx = round((nPulseCancel+1)/2);
mti_filter(midIdx) = mti_filter(midIdx) + 1;
sif = convn(sif,mti_filter,'same'); 

% display the MTI results
figure;
AbsSifToPlot = abs(sif)/max(max(abs(sif)));
imagesc(dataRange,(1:numPulses)*Tp*2,20*log10(AbsSifToPlot)); 
ylabel('Time (s)');
xlabel('Range (m)');
title('RTI with MTI clutter rejection');
colormap(jet(256));
colorbar;
axis xy; 

% ISAR imaging

z = sif;

PRF = (1/(2*Tp));disp([ 'PRF = ' num2str(PRF) ' Hz']);%Hz

transmit_wavelength = 299e6/(2.445e9 + BW/2); disp(['Wavelength = ' num2str(transmit_wavelength) ' m']);

angular_velocity = (2*pi/60)*turn_table_speed;  disp(['Angular velocity = ' num2str(angular_velocity) ' rad/s']); %rad/s

CPI = no_pulses_to_consider*(Tp*2); disp([ 'CPI = ' num2str(CPI) ' s']);%s

small_integration_angle = CPI*angular_velocity;disp(['Integration angle = ' num2str(rad2deg(small_integration_angle)) ' deg']);%rad

% Generate multiple ISAR images                                                                        

total_num_pulses = size(z,1); disp(['Total number of pulses = ' num2str(total_num_pulses)]);                                               
num_frames = round((total_num_pulses - no_pulses_to_consider)/(no_pulses_to_consider/2))-1; disp(['Number of frames = ' num2str(num_frames)]);       

for count_frames = 1:num_frames                                      
    N = (count_frames-1)*no_pulses_to_consider/2; 
	
    small_integration_window_data = z((1+N):(no_pulses_to_consider+N),:);
    small_integration_window_data = small_integration_window_data(:,:);
    
    max_value = max(max(abs(small_integration_window_data))); 
    abs_norm_small_integration_window_data = abs(small_integration_window_data)/max_value;   
 
    % Obtain range profiles within CPI
    num_range_bins = size(small_integration_window_data,2);
    window_matrix = repmat(hamming(no_pulses_to_consider),1,num_range_bins);      
    
    range_dop_map = fftshift(fft(small_integration_window_data.*window_matrix,[],1),1);      
    range_dop_map = range_dop_map(:,:);
    
    max_value = max(max(abs(range_dop_map)));                                       
    abs_norm_range_dop_map = abs(range_dop_map)/max_value;                         
    
    cross_range_resolution = transmit_wavelength/(2*small_integration_angle);disp(['Cross Range Resolution = ' num2str(cross_range_resolution) ' m']);
    cross_range = ((-no_pulses_to_consider/2:1:(no_pulses_to_consider/2 -1))*(cross_range_resolution));
    
    doppler_freq_range = ((-no_pulses_to_consider/2:1:(no_pulses_to_consider/2 -1))*(PRF/no_pulses_to_consider));
    
	% plot doppler frequency vs range
	figure(91); imagesc(dataRange,doppler_freq_range ,20*log10(abs_norm_range_dop_map)); colormap(jet(256)); h=colorbar; ylabel(h,'SNR (dB)','FontSize',12,'Rotation',270); h.Label.Position(1) = 3; axis xy; 
    title(['ISAR image number ' num2str(count_frames) ' of ' num2str(num_frames )]); ylabel('Doppler frequency (Hz)'); xlabel('Range (m)');
    
	% plot cross range vs slant range
	figure(92); imagesc(dataRange,cross_range,20*log10(abs_norm_range_dop_map)); colormap(jet(256)); h=colorbar; ylabel(h,'SNR (dB)','FontSize',12,'Rotation',270); h.Label.Position(1) = 3; axis xy;    
    title(['ISAR image number ' num2str(count_frames) ' of ' num2str(num_frames )]); ylabel('Cross-range (m)'); xlabel('Range (m)');
    
	pause(0.2);
    
end

% ---- standard DSP helper functions below ----

function [y] = fft_interp(x,M)
% perform approximate bandlimited interpolation of x by a factor of M
L = 4;
winInds = (-L*M : L*M).'/M * pi;

% get the ideal antialiasing filter's impulse response of length 2*M + 1 
winInds(L*M + 1) = 1;
myWin = sin(winInds) ./ winInds;
myWin(L*M + 1) = 1;

% use the window method; apply a hann window
myWin = myWin .* hann_window(2*L*M + 1);

% insert zeros in data and apply antialias filter via FFT
nFFT = numel(x) * M;
if isreal(x)
    y = ifft( fft(myWin,nFFT) .* repmat(fft(x),[M 1]), 'symmetric');
else
    y = ifft( fft(myWin,nFFT) .* repmat(fft(x),[M 1]) );
end
y = y([L*M+1:end 1:L*M]);

function [w] = hann_window(N)
% create a hann (cosine squared) window
w = .5 + .5*cos(2*pi*((1:N).'/(N+1) - .5));
