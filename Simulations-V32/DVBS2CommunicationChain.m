function [BER, SER, tx_symb, rx_symb] = DVBS2CommunicationChain(mod_input, Nbps_input, EbN0_ratio_input, CFO_input, t0_input, phi0_input)
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 2 - Time and Frequency Syncrhonisation


%% DVBS2CommunicationChain
% INPUTS
%       'mod_input' - Digital Modulation [PSK, PAM, QAM, Cross-QAM, Optimal-8QAM]
%       'Nbps_input' - Number of bits per symbol
%       'EbN0_ratio_input' - SNR per bit
% OUTPUTS
%       'BER' - Bit Error Ratio
%       'SER' - Symbol Error Ratio
%       'tx_symb' - Transmitted Symbols
%       'rx_symb' - Received Symbols


%% Implemented Features
% 1. [TX]   - GENERATING INFORMATION - RANDOM BITSTREAM
% 2. [TX]   - MAPPING BITSTREAM TO SYMBOLS - SUPPORTED MODULATION SCHEMES: [PSK, PAM, QAM, Cross-QAM, Optimal-8QAM] - Works for any Nbps
% 3. [TX]   - UPSAMPLING
% 4. [TX]   - FILTERING WITH ROOT-RAISED-COSINE FILTER
% 5. [CHAN] - ADDING AWGN IN THE IDEAL CHANNEL
% 6. [RX]   - SYNCHRONISATION ERRORS
% 7. [RX]   - FILTERING WITH ROOT-RAISED-COSINE FILTER
% 8. [RX]   - TIME SHIFT ON SAMPLING
% 9. [RX]   - DOWNSAMPLING
% 10.[RX]   - DEMAPPING SYMBOLS TO BITSTREAM - MAXIMUM LIKELIHOOD CRITERION
% X. [PLOT] - RELEVANT GRAPHS AND VALUES


%% IF YOU WANT TO AUTO-RUN THIS FUNCTION           - DECOMMENTS THE 4 NEXT LINES.
%% IF YOU WANT TO RUN THIS FUNCTION FROM 'FoM'     - COMMENTS THE 4 NEXT LINES.
% clear; clc; close all;
% mod_input = 'QAM';
% Nbps_input = 4;
% EbN0_ratio_input = 5;
%% RX settings 
% 1. On/Off Gardner => 0/1
Gardner_S = 1;

%% ***** RCF PARAMETERS *****
RRCF.fcutoff = 1e6; % 3dB Cutoff Frequency of the RCF [1MHz]
RRCF.beta = 0.3; % Roll-off factor of the RCF [0.3]
RRCF.M = 40; % [!] Deduct the Upsampling Factor -> Needed in order to satisfy the ISI Nyquist Criterion. Big M allows to simulate a Time Shift on sampling
RRCF.fs = RRCF.M*2*RRCF.fcutoff; % = M*fsymb - Deduct the Sampling Frequency [Hz]
RRCF.Ntaps = 16*RRCF.M + 1; % Number of taps of the RCF

%% ***** MODULATION SCHEME PARAMETERS  *****
Modu.mod = mod_input; % Select the Digital Modulation among [PSK, PAM, QAM, Cross-QAM, Optimal-8QAM] (NB: The 'Optimal-8QAM' is designed for 8 Symbols only (=> Nbps = 3))
Modu.Nbps = Nbps_input; % Select the Number of Bits Per Symbol
Modu.fsymb = 2*RRCF.fcutoff; % Select the Symbolrate [2MHz] - Common to all Modulation Schemes and independent of Nbps
Modu.bps = Modu.fsymb*Modu.Nbps; % Deduct the Bitrate [bits/s] - Depends on the Number of Bits carried per each Symbol (Nbps)
Modu.Tsymb = 1/Modu.fsymb; % Deduct the Period of a symbol [s] (or more rigorously [s/Symb])

%% ***** SYNC PARAMETERS *****
SYNC.CFO = CFO_input*Modu.fsymb; % Carrier Frequency Offset
SYNC.t0 = t0_input; % Sample Time Shift     - Always in [0:Tsymb] (Physically) - [0:M-1] (Digital Simulation)   [Physically = Tsymb*t0/M-1 = t0*Tsample]
SYNC.phi0 = phi0_input; % Carrier Phase Error - Always in [0:2pi]


%% 1. [TX] GENERATING INFORMATION - RANDOM BITSTREAM
tx_len = 10000; % Length of the bitstream
tx_bin = randi([0 1], tx_len,1); % Bitstream - Generate binary sequence (0 and 1 are equiprobable)


%% 2. [TX] MAPPING BITSTREAM TO SYMBOLS
if mod(tx_len,Modu.Nbps) ~= 0 % ZERO PADDING -> Allows to avoid to get an error if the Bitstream length is not a multiple of Nbps
    tx_bin = [tx_bin; zeros(Modu.Nbps - mod(tx_len,Modu.Nbps),1)];
    tx_len = tx_len + Modu.Nbps - mod(tx_len,Modu.Nbps);
end
tx_symb = mapping(tx_bin, Modu.Nbps, Modu.mod); % Realize the mapping according to the desired Digital Modulation Scheme and Nbps


%% 3. [TX] UPSAMPLING
tx_symb_UP = upsample(tx_symb,RRCF.M); % Upsample
tx_symb_UP_R = real(tx_symb_UP);

%% 4. [TX] FILTERING WITH ROOT-RAISED-COSINE FILTER
RRCF.IR = RRCFDesign(RRCF.beta, RRCF.Ntaps, RRCF.fs, Modu.Tsymb); % Impulse Response of the RRCF
tx_symb_UP_F = conv(tx_symb_UP,RRCF.IR); % Apply the RRCF
tx_symb_UP_F_R = real(tx_symb_UP_F);

%% 5. [CHAN] ADDING AWGN IN THE IDEAL CHANNEL
% In practice, we send in the channel the Narrowband Bandpass Signal (the Baseband Signal modulated at a high frequency)
% In this simulation, we directly send the Baseband Signal (Complex Envelope) in the channel. This is done in order to save computing power
% However, we have to take into account the fact that the energy of a Bandpass Signal is equal to the energy of its Complex Envelope divided by 2
Es_CE = (trapz(abs(tx_symb_UP_F).^2))*(1/(RRCF.fs)); % Basedband Signal Energy (Complex Envelope of the Bandpass Signal)
Es_BP = Es_CE/2; % Bandpass Signal Energy
Eb = Es_BP/tx_len; % Average energy of a bit
EbN0_dB = EbN0_ratio_input; % Energy per Bit to Noise PSD ratio [dB]
EbN0_ratio = 10^(EbN0_dB/10); % Energy per Bit to Noise PSD ratio
N0 = Eb/EbN0_ratio; % Noise PSD - We consider AWGN, so N0 is constant on the whole frequency domain
NoisePower = 2*N0*RRCF.fs; % Noise Power
noise_AWGN = sqrt(NoisePower/2)*(randn(length(tx_symb_UP_F), 1) + 1i*randn(length(tx_symb_UP_F), 1)); % Generate Zero-Mean AWGN
rx_symb_UP_F = tx_symb_UP_F + noise_AWGN; % Add the noise to the transmitted signal


%% 6. [RX] SYNCHRONISATION ERRORS
t_axis = 0:1/RRCF.fs:(1/RRCF.fs*(length(rx_symb_UP_F)-1));
exp_async = exp(1j* (2*pi*SYNC.CFO*t_axis' + SYNC.phi0) );
rx_symb_UP_F_async = exp_async.*rx_symb_UP_F;


%% 7. [RX] FILTERING WITH ROOT-RAISED-COSINE FILTER
rx_symb_UP_async = conv(rx_symb_UP_F_async,RRCF.IR); % Apply the matched RRCF
rx_symb_UP_async = rx_symb_UP_async(RRCF.Ntaps:length(rx_symb_UP_async) - RRCF.Ntaps+1); % Removes the irrelevant values


%% 8. [RX] TIME SHIFT ON SAMPLING 
% Time shift with constant offset t_0 [after the matched filter, before the down sampler]

rx_symb_UP_async_t0 = rx_symb_UP_async(SYNC.t0+1:end,1); % Thanks to the big oversampling M, we can simulated a Time Shift On Sampling
rx_symb_UP_async_t0 = [rx_symb_UP_async_t0; zeros(SYNC.t0 - 1,1)]; % Zero-Padding to match the length

rx_symb_UP_async_t0_R = real(rx_symb_UP_async_t0);


%% 9. [RX] DOWNSAMPLING
rx_symb = downsample(rx_symb_UP_async_t0,RRCF.M);

%% 10. [RX] DEMAPPING SYMBOLS TO BITSTREAM - MAXIMUM LIKELIHOOD CRITERION
rx_bin = demapping(rx_symb,Modu.Nbps,Modu.mod); % Realize the mapping according to the desired Digital Modulation Scheme and Nbps by applying the ML Criterion


%% X. [PLOT] RELEVANT GRAPHS AND VALUES
% See the function 'FoM' for all the preset Graphs
rx_bin_len = size(rx_bin,1);
BER = 1 - sum(rx_bin(1:rx_bin_len) == tx_bin(1:rx_bin_len))/(rx_bin_len); % Bit Error Ratio
SER = 1 - sum( bi2de(fliplr(reshape(rx_bin,Modu.Nbps,rx_bin_len/Modu.Nbps)')) == bi2de(fliplr(reshape(tx_bin(1:rx_bin_len),Modu.Nbps,rx_bin_len/Modu.Nbps)'))) / (rx_bin_len/Modu.Nbps); % Symbol Error Ratio


%% Gardner algorithm -> ML estimation for t_0 and correct the rx_symb
if Gardner_S == 1
    
%[T]
figure
plot(tx_symb_UP_R(200:300),'x-'); hold on
%plot(tx_symb_UP_F_R(200:300),'o-r'); hold on;
plot(rx_symb_UP_async_t0_R(200:300),'o-.k');

[est_error, y_corr] = gardner_est(rx_symb_UP_async_t0, RRCF.M, 0.1);
figure
plot(est_error);
title(['Applying Gardner to ',num2str(t0_input),' time shift unit for', num2str(2^Nbps_input), mod_input]); 
xlabel('Time (pts)'); ylabel('Estimated Error'); hold on;

rx_symb_c = downsample(y_corr,RRCF.M);
rx_bin_c = demapping(rx_symb_c,Modu.Nbps,Modu.mod); % Realize the mapping according to the desired Digital Modulation Scheme and Nbps by applying the ML Criterion
rx_bin_c_len = size(rx_bin_c,1);
BER_c = 1 - sum(rx_bin_c(1:rx_bin_c_len) == tx_bin(1:rx_bin_c_len))/(rx_bin_c_len); % Bit Error Ratio

BER
BER_c

end

end
