% Modulation and Coding Projcect
% Part 1
% 2015/Mar/19
% TEAM:MOY Michael,Omar,Yu

%% Readme
% For part one, we implement following funcions:
% 1. generate informaiton (binary) sequence
% 2. mapping the binary seq
% 3. apply filter for tx side
% 4. add Gausion noise

%% 1.[TX] Initialization and Generating Bit-stream
clear;
close all;

% generate input binary sequence by equal possiblity 0 and 1
tx_len = 4000;
Modu.mod = 'qam'; %[!] using QAM
Modu.Nbps = 4; %[!] 16QAM(4bits)
Modu.bps = 1e7; % [bits/sec] assuming 1sec for such length of bits
Modu.fsymb = Modu.bps/Modu.Nbps; % [sample/sec] 1e4
Modu.tsymb = 1/Modu.fsymb; % [sec]T:symbol time
Modu.M = 4; %[!]Up sampling factor
Modu.fsample = Modu.M*Modu.fsymb; % define sampling time 4*fsymbol=4e4

tx_bin = randi([0 1],tx_len,1);
%[tx_symb] = mapping(tx_bin,Modu.Nbps,Modu.mod); % mapping the binary seq get Re,Im
[tx_symb] = mappinga(tx_bin,Modu.Nbps,'PSK'); % mapping the binary seq get Re,Im
tx_symb_R = real(tx_symb);
tx_symb_I = imag(tx_symb);
% plot constellations figure
figure
scatter(tx_symb_R,tx_symb_I);
title('Constellation Mapping of 16QAM');


%% 2.1[TX] Filtering with RRC
% using self-define function to design RRC filter
RRCF.beta = 0.5;
RRCF.Ntaps = 33; %[?] How to select 
RRCF.td = RRCfilterDesign(RRCF.beta,RRCF.Ntaps,Modu.tsymb,Modu.fsample);
%[RRCF.td]=rrc(Modu.fsample,Modu.fsymb,RRCF.Ntaps);

% using rcosdesign ToolBox to design the filter
span = 8;
sps = 4;
% [!]Select to use system design filter or Self-Design Filter
% h1 = rcosdesign(RRCF.beta,span,sps,'sqrt');  
h1 = RRCF.td;
figure
stem(h1,'r');
hold on
stem(RRCF.td,'b');
legend('Toolbox','SelfDefine');
title('Square Raised Cosine Filter')

%RRCF.td = RRCfilterDesign(RRCF.beta,RRCF.Ntaps,1e-6,4e6); %[?] How to do Upsampling

%% 2.2[TX] Upsampling
% Upsampling *M, is to upsample the fsymbol to fsample by factor [4]
tx_usymb = upsample(tx_symb,Modu.M);

tx_usymb_R = real(tx_usymb);
tx_usymb_I = imag(tx_usymb);

% apply the filter
tx_af = conv(tx_usymb,h1);
tx_af_R = real(tx_af);

%% 3.[Channel] Adding noise

% Calculate the energy of the signa

EbN0 = 10;
EbN0mag = 10^(EbN0/10);
SignalEnergy = (trapz(abs(tx_af).^2))*(1/(Modu.fsample));
Eb = SignalEnergy/tx_len;
Eb = Eb/2;
N0 = Eb/EbN0mag;
NoisePower = 2*N0*Modu.fsample;
noise = (sqrt(NoisePower/2)*(rand(1,length(tx_af))+1i*rand(1,length(tx_af))))';
% add noise to the signal
rx_usymb = tx_af+noise;
%rx_usymb = awgn(tx_af,SNR);

% Match filter with the RRC
rx_af = conv(rx_usymb,h1);

rx_af_trunc = rx_af(RRCF.Ntaps:length(rx_af)-RRCF.Ntaps+1);
rx_af_R_trunc = real(rx_af_trunc);

% Downsampling the signal with factor M
rx_symb = downsample(rx_af_trunc,Modu.M);

% plot the tx and rx filtered signal
figure
plot(tx_usymb_R(100:200),'ro');
hold on
plot(rx_af_R_trunc(100:200),'b-x');
title('[TX] after filter')


%% 4.[RX] Demapping from Binary Sequence after the channel tx

[rx_bin] = demapping(rx_symb,Modu.Nbps,Modu.mod);

rx_symb_R = real(rx_symb);
rx_symb_I = imag(rx_symb);

% plot constellations figure
figure
scatter(rx_symb_R,rx_symb_I);
title('Constellation Mapping of 16QAM');

rtcmp = rx_bin==tx_bin;
num_rtcmp = sum(rtcmp);
rat_rtcmp = 1-num_rtcmp/tx_len % output the error ration of the coding and modulation.