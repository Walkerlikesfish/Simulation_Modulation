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

function [est_error, y_corr]=gardner_est(rx_input, M, K)

    % Init the parameters
    rx_len = length(rx_input);
    n_sym = round(rx_len/M);
    est_error = zeros(1,(n_sym)); % initialize the error mat
    y_corr = rx_input; %init the corrected rx sequence
    
    % iterate through the length of the rx
    for ii = 1:n_sym-2
        % update the error
        est_error(ii+1) = est_error(ii) + K*real(y_corr((2*ii-1)*M/2+1)...
            *(conj(y_corr(ii*M+1))-conj(y_corr((ii-1)*M+1))));
    
        % update the corrected sequence
%         y_corr((2*ii+1)*M/2+1) = interp1([ii*M+1 (2*ii+1)*M/2+1 (ii+1)*M+1],...
%             [rx_input(ii*M+1) rx_input((2*ii+1)*M/2+1) rx_input((ii+1)*M+1)],...
%             (2*ii+1)*M/2+1-round(est_error(ii+1)));
%         y_corr((ii+1)*M+1) = interp1([(2*ii+1)*M/2+1 (2*ii+2)*M/2+1 (2*ii+3)*M/2+1],...
%             [rx_input((2*ii+1)*M/2+1) rx_input((2*ii+2)*M/2+1) rx_input((2*ii+3)*M/2+1)],...
%             (ii+1)*M+1-round(est_error(ii+1)));
        y_corr((2*ii+1)*M/2+1) = interp1([ii*M+1,ii*M+2,(2*ii+1)*M/2+1,(ii+1)*M+1],...
            [rx_input(ii*M+1),rx_input(ii*M+2), rx_input((2*ii+1)*M/2+1), rx_input((ii+1)*M+1)],...
            (2*ii+1)*M/2+1-round(est_error(ii+1)));
        y_corr((ii+1)*M+1) = interp1([(2*ii+1)*M/2+1, (2*ii+2)*M/2+1, (2*ii+3)*M/2+1],...
            [rx_input((2*ii+1)*M/2+1), rx_input((2*ii+2)*M/2+1), rx_input((2*ii+3)*M/2+1)],...
            (ii+1)*M+1-round(est_error(ii+1)));
    end

    
end