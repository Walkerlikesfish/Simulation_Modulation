function [RRCf_time] = RRCfilterDesign(beta, NRRCtaps, tsymb, fsample)
%INPUTS:
%   -beta:     roll-off coefficient
%   -NRRCtaps: Number of RRC taps
%   -Tsym:     sampling rate
%OUTPUTS:
%   -RRCf_time: filter in time domain

    fmax = (1/NRRCtaps)*fsample*(NRRCtaps-1)/2;
    f = linspace(-fmax, fmax, NRRCtaps); % build the f x-axis
    Hrcf=zeros(1,NRRCtaps);
    
    for i=1:NRRCtaps
        if abs(f(i))<=((1-beta)/(2*tsymb))
            Hrcf(i) = tsymb;
        elseif abs(f(i))>((1-beta)/(2*tsymb)) && abs(f(i))<=((1+beta)/(2*tsymb))
            Hrcf(i) = tsymb/2*(1+cos(pi*tsymb/beta*(abs(f(i))-(1-beta)/(2*tsymb))));
        else 
            Hrcf(i) = 0;
        end
    end
      
    %Hrcf = fftshift(Hrcf);
    Hrrcf = sqrt(Hrcf);
    Hrrct = fftshift(ifft(ifftshift(Hrrcf), 'symmetric')); 
    normCoeff = sqrt(max(abs(conv(Hrrct,Hrrct))));
    Hrrct=Hrrct/normCoeff;
    %Hrrct = Hrrctv;
%     figure
%     stem(Hrrct);
%     title('Hrrct')
    
    %Hrctv = ifftshift(ifft(Hrcf,'symmetric'))
    %Hrct = circshift(Hrctv',-1);
    %fa = 1/Hrct((NRRCtaps+1)/2);
    %Hrct = Hrct.*fa;
    %Hrrct = Hrct.^(-2);
    
%     
%     Hrrcf = sqrt(Hrcf);
%     Hrrcf = fftshift(Hrrcf);
% %     figure
% %     plot(Hrcf);
% %     title('H(f) Square Root Cosine Filter')
%     
%     %Hrrct = ifft(Hrrcf);
%     Hrrctv = ifftshift(ifft(Hrrcf,'symmetric'));
%     Hrrct = circshift(Hrrctv',-1);
%     
    
%     figure
%     stem(Hrrct);
%     title('H(t)');
     RRCf_time=Hrrct;

end