function [rx_bin] = demapping(rx_symb,Nbps,modulation)
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 1 - Optimal communication chain over the ideal channel


%% DEMAPPING
% INPUTS
%       'rx_symb' - Symbols
%       'modulation' - Digital Modulation [PSK, PAM, QAM, Cross-QAM, Optimal-8QAM]
%       'Nbps' - Number of bits per symbol
% OUTPUT
%       'rx_bit' - Bitstream


%% INITIALIZATION
Nsymb = size(rx_symb,1); % Number of symbols

switch modulation
        %% ***** PULSE AMPLITUDE MODULATION *****
    case 'PAM'
        % Symbol to Gray-Integer
        rx_int =  maxLikelihood(rx_symb,Nbps,'PAM'); % Apply the Maximum Likelihood Criterion in order to find the exact Symbols
        % Gray-Integer to Integer
        rx_int = gray2bin(rx_int,'PAM',2^Nbps);
        % Integer to Binary
        rx_bin  = fliplr(de2bi(rx_int));
        rx_bin = reshape(rx_bin',Nsymb*Nbps,1);
        
        
        %% ***** QUADRATURE AMPLITUDE MODULATION *****
    case 'QAM'
        % Symbol to Gray-Integer
        rx_int =  maxLikelihood(rx_symb,Nbps,'QAM'); % Apply the Maximum Likelihood Criterion in order to find the exact Symbols
        % Gray-Integer to Integer
        switch Nbps
            case 3
                gray2bin_matrix = [0 1 3 2 4 5 7 6];
                for index = 1:length(rx_int)
                    rx_int(index) = gray2bin_matrix(rx_int(index)+1);
                end
            case 5
                gray2bin_matrix = [0 1 3 2 7 6 4 5 8 9 11 10 15 14 12 13 24 25 27 26 31 30 28 29 16 17 19 18 23 22 20 21];   
                for index = 1:length(rx_int)
                    rx_int(index) = gray2bin_matrix(rx_int(index)+1);
                end
            otherwise
                rx_int = gray2bin(rx_int,'QAM',2^Nbps);
        end
        % Integer to Binary
        rx_bin  = fliplr(de2bi(rx_int));
        rx_bin = reshape(rx_bin',Nsymb*Nbps,1);
        
        
        %% ***** PHASE SHIFT KEYING *****
    case 'PSK'
        % Symbol to Gray-Integer
        rx_int =  maxLikelihood(rx_symb,Nbps,'PSK'); % Apply the Maximum Likelihood Criterion in order to find the exact Symbols
        % Gray-Integer to Integer
        rx_int = gray2bin(rx_int,'PSK',2^Nbps);
        % Integer to Binary
        rx_bin  = fliplr(de2bi(rx_int));
        rx_bin = reshape(rx_bin',Nsymb*Nbps,1);
        
        
        %% ***** CROSS-QAM *****
        % In particular, considering the 'Cross-QAM', we choose to implement the '4 symbols per orbit' also called 'Star-QAM'(one of the most famous flavor of the Cross-QAM)
    case 'Cross-QAM'
        % Symbol to Gray-Integer
        rx_int =  maxLikelihood(rx_symb,Nbps,'Cross-QAM'); % Apply the Maximum Likelihood Criterion in order to find the exact Symbols
        % Gray-Integer to Integer
        rx_int = gray2bin(rx_int,'PSK',2^Nbps);
        % Integer to Binary
        rx_bin  = fliplr(de2bi(rx_int));
        rx_bin = reshape(rx_bin',Nsymb*Nbps,1);
        
        
        %% ***** OPTIMAL-8QAM *****
    case 'Optimal-8QAM'
        % Symbol to Gray-Integer
        rx_int =  maxLikelihood(rx_symb,Nbps,'Optimal-8QAM'); % Apply the Maximum Likelihood Criterion in order to find the exact Symbols
        % Gray-Integer to Integer
        rx_int = gray2bin(rx_int,'PSK',2^Nbps);
        % Integer to Binary
        rx_bin  = fliplr(de2bi(rx_int));
        rx_bin = reshape(rx_bin',Nsymb*Nbps,1);
end
end