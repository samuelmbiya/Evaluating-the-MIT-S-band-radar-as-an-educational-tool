%% Processing the recordings for ISAR imaging

clear all;
close all;

% Radar parameters

CPI = 0.8;%s
Range_profiles_considered = 16;
PRF = 20;%Hz
wavelength = 3e8/2.4e9;% m

% Turntable parameters

Rotational_rate = 0.380;% rad/s change to 0.676 for fast rotation speed
R = 0.5;
Velocity = R*Rotational_rate;% m/s
Maximum_velocity = 2*Velocity/wavelength;

% ISAR image resolution

Doppler_resolution = 1/CPI;% Hz
Cross_range_resolution = wavelength/(2*Rotational_rate*CPI)% m

% ISAR recording dataset
All_ISAR_wav_recordings = {'1_corner_reflectors_1m_downrange_slow_FMCW_test2.wav';'2_corner_reflectors_1m_downrange_slow.wav'};

% Read in the recording
ISAR_recording_no = 1;

mit_radar_isar_processing(All_ISAR_wav_recordings{ISAR_recording_no}); 
