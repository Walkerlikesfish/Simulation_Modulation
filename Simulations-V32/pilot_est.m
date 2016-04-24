% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 2 - Time and Frequency Syncrhonisation
% Gardner Function
% carry out Gardner to Maximum Likelyhood estimation of t_0
% input:
%   -sig_in: receiver rx_input downsampled
%   -pilot: pilot symbols
%   -K: constant in the Gardner algorithm
% output:
%   -n_est: estimated best synchornized time
%   -CFO_est: estimated CFO

function [n_est,CFO_est] = pilot_est(sig_in,pilot,f_sym)

N = size(pilot,1);     % N is the size of the Pilot symbols
K = 16;                % Cross-Correlation Window [!] Select to improve performance
Tsymb = 1/f_sym;       % Symbol interval
len = length(sig_in);  % length of the whole input signal

matD = zeros(K,len-N+1); % initialize the diff corr matrix

for kk=1:K   % iterating the corr diff
    for l=kk:N-1 % going through the points
        m = (conj(sig_in(l:len-N+l))*pilot(l)...
            .* (conj(conj(sig_in(1+l-kk:len-N+l-kk+1))*pilot(l-kk+1))));
        matD(kk,:)= matD(kk,:) + m.';
    end
    matD(kk,:)=matD(kk,:)/(N-kk);
end

tot_matD=sum(abs(matD));
[~,n_est]= max(tot_matD);

CFO_est=0;

for kk=1:K
    CFO_est= CFO_est + phase(matD(kk,n_est))/(2*pi*kk*Tsymb);
end

CFO_est=-CFO_est/K; % the estimation 

end