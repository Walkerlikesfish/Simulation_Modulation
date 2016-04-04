% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 1 - Optimal communication chain over the ideal channel


%% Implemented Features
% For Part 1, we implemented the following functions
% 1. [TX] GENERATE INFORMATION - RANDOM BITSTREAM
% 2. [TX] MAP BITSTREAM TO SYMBOLS
% 3. Apply filter for tx side
% 4. Add AWGN
clear; clc; close all;


%% ***** MODULATOR PARAMETERS  ***** 
Modu.mod = 'QAM'; % Select the Digital Modulation among [PSK, PAM, QAM, Cross-QAM]
Modu.Nbps = 6; % Select the number of bits per symbol among [1, 2, 3, 4, 5, 6]
Modu.Ns = 2^Modu.Nbps; % Number of symbols = 2^Nbps [2, 4, 8, 16, 32, 64]
Modu.bps = 1e7; % Select the Bitrate [bits/s]
Modu.fsymb = Modu.bps/Modu.Nbps; % Symbolrate [Symb/s]
Modu.tsymb = 1/Modu.fsymb; % Period of a symbol [s/Symb]
Modu.M = 4; % Select the Upsampling Factor -> Allows to emulate a higher rate sampling -> Needed in order to satisfy the ISI Nyquist Criterion
Modu.fsample = Modu.M*Modu.fsymb; % Sampling Frequency [Hz] (or more rigorously [Symb/s])


%% 1. [TX] GENERATE INFORMATION - RANDOM BITSTREAM
tx_len = 2000; % Length of the bitstream
tx_bin = randi([0 1],tx_len,1); % Bitstream - Generate binary sequence (0 and 1 are equiprobable)

if mod(tx_len,Modu.Nbps) ~= 0 % ZERO PADDING -> Allows to avoid to get an error if the Bitstream length is not a multiple of Nbps
    tx_bin = [tx_bin; zeros(Modu.Nbps - mod(tx_len,Modu.Nbps),1)];
    tx_len = tx_len + Modu.Nbps - mod(tx_len,Modu.Nbps);
end

%% 2. [TX] MAP BITSTREAM TO SYMBOLS
[tx_symb] = mapping(tx_bin,Modu.Nbps,Modu.mod); % Realize the mapping - Digital Modulation
tx_symb_R = real(tx_symb);
tx_symb_I = imag(tx_symb);

% Plot the Constellation Symbols of the selected mapping
fig1 = figure; scatter(tx_symb_R,tx_symb_I);
ax=gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; % Set the axis to the origin
ax.XLim = [min(tx_symb_R) - 0.25, max(tx_symb_R) + 0.25]; ax.YLim = [min(tx_symb_R) - 0.25, max(tx_symb_R) + 0.25]; ax.XLimMode = 'manual'; axis square; % Set the axis limits
title(['Constellation Symbols of ', num2str(Modu.Ns), Modu.mod]); xlabel('Real'); ylabel('Imaginary'); hold on;

%% 3.1 [TX] Filtering with RRC
% using self-define function to design RRC filter
RRCF.beta = 0.5; % Set the value of the Roll-off factor (caracterises the steepness of a transmission function with frequency)
RRCF.Ntaps = 33; %[?] How to select 
RRCF.td = RRCfilterDesign(RRCF.beta,RRCF.Ntaps,Modu.tsymb,Modu.fsample);
%[RRCF.td]=rrc(Modu.fsample,Modu.fsymb,RRCF.Ntaps);
% using rcosdesign ToolBox to design the filter
span = 8;
sps = 4;
% h1 = rcosdesign(RRCF.beta,span,sps,'sqrt');
h1 = RRCF.td;

% Plot the Square Raised Cosine Filter
fig2 = figure; stem(h1,'r'); 
title(['Square Raised Cosine Filter (Tsymb = ', num2str(Modu.tsymb*10^6), 'µs/symb and fs = ', num2str(Modu.fsample/(10^6)), 'MHz)']); xlabel('Samples'); ylabel('Amplitude'); legend('Toolbox'); hold on;

%RRCF.td = RRCfilterDesign(RRCF.beta,RRCF.Ntaps,1e-6,4e6); %[?] How to do Upsampling

%% 3.2[TX] Upsampling
% Upsampling *M, is to upsample the fsymbol to fsample by factor [4]
tx_usymb = upsample(tx_symb,Modu.M);

tx_usymb_R = real(tx_usymb);
tx_usymb_I = imag(tx_usymb);

% apply the filter
tx_usymb_af = conv(tx_usymb,h1);
tx_usymb_af_R = real(tx_usymb_af);
tx_usymb_af_I = imag(tx_usymb_af);

%% 4.[Channel] Adding noise

% Calculate the energy of the signal

EbN0 = 10;
EbN0mag = 10^(EbN0/10);
SignalEnergy = (trapz(abs(tx_usymb_af).^2))*(1/(Modu.fsample));
Eb = SignalEnergy/tx_len;
Eb = Eb/2;
N0 = Eb/EbN0mag;
NoisePower = 2*N0*Modu.fsample;
noise = rand(1,length(tx_usymb_af)) .* sign(rand(length(tx_usymb_af),1)-0.5)' + 1i*rand(1,length(tx_usymb_af)) .* sign(rand(length(tx_usymb_af),1)-0.5)';
noise = (sqrt(NoisePower/2)*noise)';

% add noise to the signal
rx_usymb = tx_usymb_af + noise;
rx_usymb_R = real(rx_usymb);
rx_usymb_I = imag(rx_usymb);
%rx_usymb = awgn(tx_af,SNR);

% Match filter with the RRC
rx_af = conv(rx_usymb,h1);
rx_usymb_af_R = real(rx_af);
rx_usymb_af_I = imag(rx_af);

rx_af_trunc = rx_af(RRCF.Ntaps:length(rx_af)-RRCF.Ntaps+1);
rx_af_R_trunc = real(rx_af_trunc);
rx_af_I_trunc = imag(rx_af_trunc);

% Downsampling the signal with factor M
rx_symb = downsample(rx_af_trunc,Modu.M);

% Plot the ... %% TODO modify this part !
fig3 = figure; 
subplot(2,2,1); plot(tx_usymb_R(100:200),'ro'); hold on; plot(rx_usymb_R(100:200),'b-x');
title('Real Part of the Signals (before RRC Filter)'); xlabel('Samples'); ylabel('Amplitude'); legend('Real Part of the Transmitted Signal', 'Real Part of the Received Signal (before RRC Filter)'); hold on;
subplot(2,2,2); plot(tx_usymb_I(100:200),'ro'); hold on; plot(rx_usymb_I(100:200),'b-x');
title('Imaginary Part of the Signals (before RRC Filter)'); xlabel('Samples'); ylabel('Amplitude'); legend('Imaginary Part of the Transmitted Signal', 'Imaginary Part of the Received Signal (before RRC Filter)'); hold on;
subplot(2,2,3); plot(tx_usymb_af_R(100:200),'ro'); hold on; plot(rx_usymb_af_R(100:200),'b-x');
title('Real Part of the Signals (after RRC Filter)'); xlabel('Samples'); ylabel('Amplitude'); legend('Imaginary Part of the Received Signal', 'Imaginary Part of the Received Signal (after RRC Filter)'); hold on;
subplot(2,2,4); plot(tx_usymb_af_I(100:200),'ro'); hold on; plot(rx_usymb_af_I(100:200),'b-x');
title('Imaginary Part of the Signals (after RRC Filter)'); xlabel('Samples'); ylabel('Amplitude'); legend('Real Part of the Received Signal', 'Real Part of the Received Signal (after RRC Filter)'); hold on;

%% 5.[RX] Demapping from Binary Sequence after the channel tx
[rx_bin] = demapping(rx_symb,Modu.Nbps,Modu.mod);
rx_symb_R = real(rx_symb);
rx_symb_I = imag(rx_symb);

% Plot the Constellation Symbols at the receiver
fig4 = figure; scatter(rx_symb_R,rx_symb_I);
ax=gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin';
ax.XLim = [min(tx_symb_R) - 0.25, max(tx_symb_R) + 0.25]; ax.YLim = [min(tx_symb_R) - 0.25, max(tx_symb_R) + 0.25]; ax.XLimMode = 'manual'; axis square;
title(['Constellation Symbols of ', num2str(Modu.Ns), Modu.mod]); xlabel('Real'); ylabel('Imaginary'); hold on;

%% BIT ERROR RATIO
num_rtcmp = sum(rx_bin == tx_bin);
BER = 1 - num_rtcmp/tx_len % Bit Error Ratio of the coding and modulation.