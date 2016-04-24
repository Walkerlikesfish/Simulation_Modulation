function [IR_RRC, H_RRCF, t_axis, f_axis] = RRCFDesign(Beta, Ntaps, fs, Tsymb)
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 2 - Time and Frequency Syncrhonisation


%% RRCFDesign
% INPUTsymb
%       'Beta' - Roll-off Factor
%       'Ntaps' - Number of RRC Taps [Needs to be odd] - Set the Time Extension and the Frequency Resolution
%       'fs' - Sampling Frequency [Hz]
%       'Tsymb' - Symbol Period [s]
% OUTPUTsymb
%       'IR_RRC' - Impulse Response of the designed RRCF
%       'H_RRCF' - Spectrum of the designed RRCF


%% The RCF is designed in the frequency domain
%% The Impulse Response is then obtained by IFFT
fmax = (1/Ntaps)*fs*(Ntaps-1)/2; % Maximal frequency on the axis - Tends to fs/2 for Ntaps tending to infinity
f_axis = linspace(-fmax, fmax, Ntaps);
t_axis = (-(Ntaps-1)/2:(Ntaps-1)/2)./(2*fmax);
H_RCF = zeros(1,Ntaps); % Fourier Transform of the RCF

for i=1:Ntaps % Ntaps seTsymb the Time Extension and the Frequency Resolution
    if abs(f_axis(i))<=((1-Beta)/(2*Tsymb))
        H_RCF(i) = Tsymb;
    elseif abs(f_axis(i))>((1-Beta)/(2*Tsymb)) && abs(f_axis(i))<=((1+Beta)/(2*Tsymb))
        H_RCF(i) = Tsymb/2*(1+cos(pi*Tsymb/Beta*(abs(f_axis(i))-(1-Beta)/(2*Tsymb))));
    else
        H_RCF(i) = 0;
    end
end

H_RRCF = sqrt(H_RCF);
IR_RRC = fftshift(ifft(ifftshift(H_RRCF), 'symmetric')); 
normCoeff = sqrt(max(abs(conv(IR_RRC,IR_RRC))));
IR_RRC=IR_RRC/normCoeff; % Normalize -> the RCF will then be equal to 1 @ t=0 (and not the RRCF !)
end