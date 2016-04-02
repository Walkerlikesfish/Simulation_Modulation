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
tx_len = 40000;
Modu.mod = 'qam'; %[!] using QAM
Modu.Nbps = 4; %[!] 16QAM(4bits)
Modu.bps = tx_len; % [bits/sec] assuming 1sec for such length of bits
Modu.fsymb = tx_len/Modu.Nbps; % [sample/sec] 1e4
Modu.tsymb = 1/Modu.fsymb; % [sec]T:symbol time
Modu.M = 4; %[!]Up sampling factor
Modu.fsample = Modu.M*Modu.fsymb; % define sampling time 4*fsymbol=4e4

tx_bin = randi([0 1],tx_len,1);
[tx_symb] = mapping(tx_bin,Modu.Nbps,Modu.mod); % mapping the binary seq get Re,Im
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

% using rcosdesign ToolBox to design the filter
span = 8;
sps = 4;
h1 = rcosdesign(RRCF.beta,span,sps,'sqrt');
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
% trunc the result
tx_af_trunc = tx_af((RRCF.Ntaps+1)/2:length(tx_af)-(RRCF.Ntaps+1)/2+1);
tx_af_R_trunc = real(tx_af_trunc);


%% 3.[Channel] Adding noise

% Calculate the energy of the signal
SignalEnergy = (trapz(abs(tx_af_R_trunc).^2))*(1/Modu.fsample);
Eb = SignalEnergy/Modu.Nbps;
Eb = Eb/2; %[?] why we need to divide by 2 here [?]
EbN0 = 3; %[!]EbN0: signal noise ratio[dB]
N0 = Eb/db2mag(EbN0);
NoisePower = 0;
noise = zeros(1,length(tx_af_trunc));
noise = sqrt(NoisePower/2)*(randn(1,length(noise))+1i*randn(1,length(noise)));

% add noise to the signal
rx_usymb = tx_af_trunc+noise';
% rx_usymb_ = tx_af;

% Match filter with the RRC
rx_af = conv(rx_usymb,h1);
% rx_af_ = conv(rx_usymb_,h1);

rx_af_trunc = rx_af((RRCF.Ntaps+1)/2:length(rx_af)-(RRCF.Ntaps+1)/2+1);
rx_af_R_trunc = real(rx_af_trunc);

% rx_af_trunc_ = rx_af_(RRCF.Ntaps:length(rx_af_)-RRCF.Ntaps+1);

% Downsampling the signal with factor M
rx_symb = downsample(rx_af_trunc,Modu.M);
% rx_symb_ = downsample(rx_af_trunc_,Modu.M);

% plot the tx and rx filtered signal
figure
plot(tx_usymb_R(1:100),'rx');
hold on
plot(tx_af_R_trunc(1:100),'g-o');
hold on
plot(rx_af_R_trunc(1:100),'b-x');
title('[TX] after filter')


%% 4.[RX] Demapping from Binary Sequence after the channel tx

[rx_bin] = demapping(rx_symb,Modu.Nbps,Modu.mod);
% [rx_bin_] = demapping(rx_symb_,Modu.Nbps,Modu.mod);

rx_symb_R = real(rx_symb);
rx_symb_I = imag(rx_symb);

% rx_symb_R_ = real(rx_symb_);
% rx_symb_I_ = imag(rx_symb_);

% plot constellations figure
figure
scatter(rx_symb_R,rx_symb_I);
title('Constellation Mapping of 16QAM');

rtcmp = rx_bin==tx_bin;
num_rtcmp = sum(rtcmp);
rat_rtcmp = num_rtcmp/tx_len % output the error ration of the coding and modulation.
