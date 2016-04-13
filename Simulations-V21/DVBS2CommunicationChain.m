function [BER, SER, tx_symb, rx_symb] = DVBS2CommunicationChain(mod_input, Nbps_input, EbN0_ratio_input)
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
% Nbps_input = 6;
% EbN0_ratio_input = 5;


%% ***** MODULATION SCHEME PARAMETERS  *****
Modu.mod = mod_input; % Select the Digital Modulation among [PSK, PAM, QAM, Cross-QAM, Optimal-8QAM] (NB: The 'Optimal-8QAM' is designed for 8 Symbols only (=> Nbps = 3))
Modu.Nbps = Nbps_input; % Select the Number of Bits Per Symbol
Modu.fsymb = 10e6; % Select the Symbolrate [Symb/s] - Common to all Modulation Schemes and independent of Nbps
Modu.bps = Modu.fsymb*Modu.Nbps; % Deduct the Bitrate [bits/s] - Depends on the Number of Bits carried per each Symbol (Nbps)
Modu.Tsymb = 1/Modu.fsymb; % Deduct the Period of a symbol [s] (or more rigorously [s/Symb])


%% ***** RCF PARAMETERS *****
RRCF.fcutoff = 1e6; % Cutoff frequency of the RCF [1MHz]
RRCF.beta = 0.3; % Roll-off factor of the RCF [0.3]
RRCF.Ntaps = 33; % Number of taps of the RCF
RRCF.M = 20; % Deduct the Upsampling Factor -> Needed in order to satisfy the ISI Nyquist Criterion. Big M allows to simulate a Time Shift on sampling
RRCF.fs = RRCF.M*Modu.fsymb; % Deduct the Sampling Frequency [Hz]


%% ***** SYNC PARAMETERS *****
SYNC.CFO = 0; % Carrier Frequency Offset
SYNC.t0 = 0; % Sample Time Shift     - Always in [0:Tsymb] (Physically) - [0:M-1] (Digital Simulation)
SYNC.phi0 = pi/4; % Carrier Phase Error - Always in [0:2pi]



%% 1. [TX] GENERATING INFORMATION - RANDOM BITSTREAM
tx_len = 5000; % Length of the bitstream
tx_bin = randi([0 1], tx_len,1); % Bitstream - Generate binary sequence (0 and 1 are equiprobable)


%% 2. [TX] MAPPING BITSTREAM TO SYMBOLS
if mod(tx_len,Modu.Nbps) ~= 0 % ZERO PADDING -> Allows to avoid to get an error if the Bitstream length is not a multiple of Nbps
    tx_bin = [tx_bin; zeros(Modu.Nbps - mod(tx_len,Modu.Nbps),1)];
    tx_len = tx_len + Modu.Nbps - mod(tx_len,Modu.Nbps);
end
tx_symb = mapping(tx_bin, Modu.Nbps, Modu.mod); % Realize the mapping according to the desired Digital Modulation Scheme and Nbps


%% 3. [TX] UPSAMPLING
tx_symb_UP = upsample(tx_symb,RRCF.M); % Upsample


%% 4. [TX] FILTERING WITH ROOT-RAISED-COSINE FILTER
RRCF.IR = RRCFDesign(RRCF.beta, RRCF.Ntaps, RRCF.fs, Modu.Tsymb); % Impulse Response of the RRCF
[RRCF.IR] = RRCfilterDesign(RRCF.beta, RRCF.Ntaps, RRCF.fs, Modu.Tsymb);
%INPUTS:
tx_symb_UP_F = conv(tx_symb_UP,RRCF.IR); % Apply the RRCF


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
%rx_symb_UP_F = tx_symb_UP_F;


%% 6. [RX] SYNCHRONISATION ERRORS
t_axis = 0:(Modu.Tsymb/(RRCF.M)):(Modu.Tsymb*(length(rx_symb_UP_F)-1)/(RRCF.M));
exp_async = exp(1j* (2*pi*SYNC.CFO*t_axis' + SYNC.phi0) );
rx_symb_UP_F_async = exp_async.*rx_symb_UP_F;


%% 7. [RX] FILTERING WITH ROOT-RAISED-COSINE FILTER
rx_symb_UP_async = conv(rx_symb_UP_F_async,RRCF.IR); % Apply the matched RRCF
rx_symb_UP_async = rx_symb_UP_async(RRCF.Ntaps:length(rx_symb_UP_async) - RRCF.Ntaps+1); % Removes the irrelevant values


%% 8. [RX] TIME SHIFT ON SAMPLING
rx_symb_UP_async_t0 = rx_symb_UP_async(SYNC.t0+1:end,1); % Thanks to the big oversampling M, we can simulated a Time Shift On Sampling
rx_symb_UP_async_t0 = [rx_symb_UP_async_t0 zeros(1, SYNC.t0 - 1)]; % Zero-Padding to match the length


%% 9. [RX] DOWNSAMPLING
rx_symb = downsample(rx_symb_UP_async_t0,RRCF.M);


%% 10. [RX] DEMAPPING SYMBOLS TO BITSTREAM - MAXIMUM LIKELIHOOD CRITERION
rx_bin = demapping(rx_symb,Modu.Nbps,Modu.mod); % Realize the mapping according to the desired Digital Modulation Scheme and Nbps by applying the ML Criterion


%% X. [PLOT] RELEVANT GRAPHS AND VALUES
% See the function 'FoM' for all the preset Graphs
BER = 1 - sum(rx_bin == tx_bin)/tx_len; % Bit Error Ratio
SER = 1 - sum( bi2de(fliplr(reshape(rx_bin,Modu.Nbps,tx_len/Modu.Nbps)')) == bi2de(fliplr(reshape(tx_bin,Modu.Nbps,tx_len/Modu.Nbps)'))) / (tx_len/Modu.Nbps); % Symbol Error Ratio
end
