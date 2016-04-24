% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 2 - Time and Frequency Syncrhonisation
% Gardner Function
% carry out Gardner to Maximum Likelyhood estimation of t_0
% input:
%   -rx_input: receiver rx_input upsampled by factor M
%   -M: upsampling factor
%   -K: constant in the Gardner algorithm
% output:
%   -est_error: estimated est_error
%   -y_corr: corrected sequence

function [CFO_esti] = pilot_est( signal, a , fsym )

N=numel(a);
K=16;
T=1/fsym;
L=length(signal);

D=zeros(K,L-N+1);
for kk=1:K
    for ll=kk:N-1
        D(kk,:)= D(kk,:)+ transpose(conj(signal(ll:L-N+ll))*a(ll) .* (conj(conj(signal(1+ll-kk:L-N+ll-kk+1))*a(ll-kk+1))));
    end
    D(kk,:)=D(kk,:)/(N-kk);
end

Dtot=sum(abs(D));
[maximum,n_est]= max(Dtot);

CFO_esti=0;
for kk=1:K
    CFO_esti= CFO_esti + phase(D(kk,n_est))/(2*pi*kk*T);
end
CFO_esti=-CFO_esti/K;
end