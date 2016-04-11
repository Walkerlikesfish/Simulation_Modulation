function [tx_symb] = mapping(tx_bin,Nbps,modulation) % Maps bitstream to Symbols depending on the selected mapping
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 1 - Optimal communication chain over the ideal channel


%% MAPPING
% INPUTS
%       'tx_bin' - Bitstream
%       'modulation' - Digital Modulation [PSK, PAM, QAM, Cross-QAM, Optimal-8QAM]
%       'Nbps' - Number of bits per symbol
% OUTPUT
%       'tx_symb' - Symbols


%% INITIALIZATION AND VALIDITY CHECK
if mod(Nbps, 1) ~=0 || Nbps <= 0
    error('You have entered an unsupported Number of Bits Per Symbol. Nbps must be a positive real integer.');
end
Nsymb = size(tx_bin,1)/Nbps; % Number of symbols
tx_bin = reshape(tx_bin,Nbps,Nsymb)'; % Reshapes the vector into a matrix: each row corresponds to a symbol (Nbps bits)

switch modulation
        %% ***** PULSE AMPLITUDE MODULATION *****
    case 'PAM'
        % Binary to Integer - Convert 'Binary Word' to corresponding 'Integer' [10 -> 2]
        tx_int = bi2de(fliplr(tx_bin));  %[!]fliplr to adjust MSB and LSB in the right position for f_bi2de
        % Integer to Gray-Integer - Convert 'Integer' to corresponding 'Gray-Integer' according to the desired Modulation Scheme and Npbs [2 <-> 3]
        tx_int = bin2gray(tx_int,'PAM',2^Nbps);
        % Gray-Integer to Symbol - Convert 'Gray Integer' to corresponding 'Symbol' according to the desired Modulation Scheme and Npbs
        distance = 2/(2^Nbps-1);
        tx_symb = distance * (tx_int - (2^Nbps-1)/2);
        
        
        %% ***** QUADRATURE AMPLITUDE MODULATION *****
    case 'QAM'
        % Binary to Integer
        tx_int = bi2de(fliplr(tx_bin));
        % Integer to Gray-Integer
        switch Nbps
            case 3
                bin2gray_matrix = [0 1 3 2 4 5 7 6];
                for index = 1:length(tx_int)
                    tx_int(index) = bin2gray_matrix(tx_int(index)+1);
                end
            case 5
                bin2gray_matrix = [0 1 3 2 6 7 5 4 8 9 11 10 14 15 13 12 24 25 27 26 30 31 29 28 16 17 19 18 22 23 21 20];
                for index = 1:length(tx_int)
                    tx_int(index) = bin2gray_matrix(tx_int(index)+1);
                end
            otherwise
                tx_int = bin2gray(tx_int,'QAM',2^Nbps);
        end
        
        % Gray-Integer to Symbol
        Nbps_I = ceil(Nbps/2); % Number of symbols carried in the 'Inphase Carrier'
        Nbps_Q = Nbps - Nbps_I; % Number of symbols carried in the 'Quadrature Carrier'
        distance = 2/(2^Nbps_I-1);
        tx_symb = distance*( (mod(tx_int,2^(Nbps_I)) - (2^Nbps_I-1)/2)  -  1i*(floor(tx_int/(2^Nbps_I)) - (2^Nbps_Q-1)/2) );
        
        
        %% ***** PHASE SHIFT KEYING *****
    case 'PSK'
        % Binary to Integer
        tx_int = bi2de(fliplr(tx_bin));
        % Integer to Gray-Integer
        tx_int = bin2gray(tx_int,'PSK',2^Nbps);
        % Gray-Integer to Symbol
        tx_symb = exp(1i*pi*(0+ 2/(2^Nbps)*tx_int));  %[?] 3/4 <- 0
        
        
        %% ***** CROSS-QAM *****
        % In particular, considering the 'Cross-QAM', we choose to implement the '4 symbols per orbit' also called 'Star-QAM' (one of the most famous flavour of the Cross-QAM)
    case 'Cross-QAM'
        switch Nbps
            case 1 % Simple PSK
                tx_int = bi2de(fliplr(tx_bin));
                tx_symb = exp(1i*pi*(3/4 + 2/(2^Nbps)*tx_int));
            otherwise
                % Binary to Integer
                tx_int = bi2de(fliplr(tx_bin));
                % Integer to Gray-Integer
                tx_int = bin2gray(tx_int,'PSK',2^Nbps);
                % Gray-Integer to Symbol
                tx_symb = zeros(length(tx_int),1);
                for index = 1:length(tx_int)
                    radius = (2 + floor(tx_int(index)/4)) / (2^Nbps/4 + 1);
                    phase = (mod(floor(tx_int(index)/4), 2) == 0)*3/4;
                    tx_symb(index) = radius*exp(1i*pi*(phase + 1/2*mod(tx_int(index),4)));
                end
        end
        
        
        %% ***** OPTIMAL-8QAM *****
    case 'Optimal-8QAM'
        if Nbps ~= 3
            error('This particular Modulation Scheme only works with 3 Bits per Symbol. For Nbps other than 3, please choose among [PSK, PAM, QAM, Cross-QAM].');
        end
        % Binary to Integer
        tx_int = bi2de(fliplr(tx_bin));
        % Integer to Gray-Integer
        tx_int = bin2gray(tx_int,'PSK',2^Nbps);
        % Gray-Integer to Symbol
        tx_symb = zeros(length(tx_int),1);
        for index = 1:length(tx_int)
            if tx_int(index) > 3
                radius = 1;
            else
                radius = sqrt(2) / (1 + sqrt(3));
            end
            phase = (floor(tx_int(index)/4) == 0)*3/4;
            tx_symb(index) = radius*exp(1i*pi*(phase + 1/2*mod(tx_int(index),4)));
        end
        
    otherwise
        error('You have entered an unsupported Digital Modulation scheme. Please choose among [PSK, PAM, QAM, Cross-QAM, Optimal-8QAM].');
end