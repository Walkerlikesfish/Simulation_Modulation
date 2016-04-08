function [tx_symb] = mapping(tx_bin,Nbps,modulation) % Maps bitstream to Symbols depending on the selected mapping
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 1 - Optimal communication chain over the ideal channel


%% MAPPING
% INPUTS
    % 'tx_bin' - Bitstream
    % 'modulation' - Digital Modulation [PSK, PAM, QAM, CQAM]
    % 'Nbps' - Number of bits per symbol [1, 2, 3, 4, 5, 6]
% OUTPUT
    % 'tx_symb' - Symbols

    
%% INITIALIZATION AND CHECK VALIDITY
if Nbps > 6 || Nbps <= 0
    error('You have entered an unsupported Number of Bits Per Symbol. Please choose among [1, 2, 3, 4, 5, 6].');
end
Nsymb = size(tx_bin,1)/Nbps; % Number of symbols
tx_bin = reshape(tx_bin,Nbps,Nsymb)'; % Reshapes the vector into a matrix: each row corresponds to a symbol (Nbps bits)

switch modulation
    %% ***** PULSE AMPLITUDE MODULATION *****
    case 'PAM'
        % Gray to binary
        for ii = 2:Nbps,
           tx_bin(:,ii) = xor( tx_bin(:,ii-1) , tx_bin(:,ii) ); 
        end
        % Binary to integer: traduction of the binary words of tx_bin into the corresponding symbols of tx_int
        tx_int = bi2de(fliplr(tx_bin)); % e.g: if Nbps = 4 -> M = 16 -> each row of tx_int can range from 0 to 15
        % Integer to symbol
        distance = 2/(2^Nbps-1); % Euclidian Distance between adjacent symbols: for PAM it's the same for all symbols couple
        tx_symb = distance * (tx_int - (2^Nbps-1)/2);

     
        
    %% ***** QUADRATURE AMPLITUDE MODULATION *****
    case 'QAM'
        switch Nbps
            case 1 % Simple PAM
                tx_int = bi2de(fliplr(tx_bin));
                distance = 2/(2^Nbps-1);
                tx_symb = distance * (tx_int - (2^Nbps-1)/2);

            otherwise % RECTANGULAR-QAM
                Nbps_I = ceil(Nbps/2); % Number of symbols carried in the 'Inphase Carrier'
                Nbps_Q = Nbps - Nbps_I; % Number of symbols carried in the 'Quadrature Carrier'

                %% *** REAL PART ***
                tx_bin_I = tx_bin(:,1:Nbps_I); % Matrix whose rows carry the first half of the bits present in each symbol (the ones appearing in the 'Inphase Carrier')
                % Gray to binary
                for ii = 2:Nbps_I,
                   tx_bin_I(:,ii) = xor( tx_bin_I(:,ii-1) , tx_bin_I(:,ii) ); 
                end
                % Binary to integer
                tx_int_I = bi2de(fliplr(tx_bin_I));
                % Integer to symbol
                distance = 2/(2^Nbps_I-1); % Euclidian Distance between adjacent symbols of the Real Part
                tx_symb_I = distance * (tx_int_I - (2^Nbps_I-1)/2);
                
                %% *** IMAGINARY PART ***
                tx_bin_Q = tx_bin(:,Nbps_I+1:end);
                % Gray to binary
                for ii = 2:Nbps_Q,
                   tx_bin_Q(:,ii) = xor( tx_bin_Q(:,ii-1) , tx_bin_Q(:,ii) ); 
                end
                % Binary to integer
                tx_int_Q = bi2de(fliplr(tx_bin_Q));
                % Integer to symbol
                tx_symb_Q = distance * (tx_int_Q - (2^Nbps_Q-1)/2);
                %% *** COMPLEX SYMBOL ***
                tx_symb = tx_symb_I + 1i*tx_symb_Q;
        end
        
        

    %% ***** PHASE SHIFT KEYING *****
    case 'PSK'
        % Gray to binary
        for ii = 2:Nbps,
           tx_bin(:,ii) = xor( tx_bin(:,ii-1) , tx_bin(:,ii) ); 
        end
        % Binary to integer
        tx_int = bi2de(fliplr(tx_bin));
        % Integer to symbol
        tx_symb = exp(1i*pi*(3/4 + 2/(2^Nbps)*tx_int));

        
        
    %% ***** CROSS-QAM ***** 
    case 'Cross-QAM' % In particular, we choose to implement the "4 symbols per orbit" Cross-QAM (one of the most famous flavor of the Cross-QAM)
        % Gray to binary
        for ii = 2:Nbps,
           tx_bin(:,ii) = xor( tx_bin(:,ii-1) , tx_bin(:,ii) ); 
        end
        % Binary to integer
        tx_int = bi2de(fliplr(tx_bin));
        tx_symb = zeros(1,length(tx_int));
        
        %% *** COMPLEX SYMBOL ***
        % Integer to symbol
        switch Nbps
            case {1,2} % One Orbit => Same case as the PSK
                tx_symb = exp(1i*pi*(3/4 + 2/(2^Nbps)*tx_int));
                
            case 3 % Two Orbits
                key_Amp =   {0, 1, 2, 3, 4, 5, 6, 7};
                value_Amp = [1/2, 1/2, 1/2, 1/2, 2/2, 2/2, 2/2, 2/2];
                Amp = containers.Map(key_Amp,value_Amp);
                key_phase =   {0, 1, 2, 3, 4, 5, 6, 7};
                value_phase = [3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0];
                Phase = containers.Map(key_phase,value_phase);
                for indexes = 1:length(tx_int)
                    tx_symb(indexes) = Amp(tx_int(indexes))*exp(1i*pi*(Phase(tx_int(indexes)) + 1/2*tx_int(indexes)));
                end
                
            case 4 % Four Orbit
                key_Amp =   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
                value_Amp = [2/5, 2/5, 2/5, 2/5, 3/5, 3/5, 3/5, 3/5, 4/5, 4/5, 4/5, 4/5, 5/5, 5/5, 5/5, 5/5];
                Amp = containers.Map(key_Amp,value_Amp);
                key_phase =   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
                value_phase = [3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0];
                Phase = containers.Map(key_phase,value_phase);
                for indexes = 1:length(tx_int)
                    tx_symb(indexes) = Amp(tx_int(indexes))*exp(1i*pi*(Phase(tx_int(indexes)) + 1/2*tx_int(indexes)));
                end
                
            case 5 % 8 Orbits
                key_Amp =   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
                value_Amp = [2/9, 2/9, 2/9, 2/9, 3/9, 3/9, 3/9, 3/9, 4/9, 4/9, 4/9, 4/9, 5/9, 5/9, 5/9, 5/9, 6/9, 6/9, 6/9, 6/9, 7/9, 7/9, 7/9, 7/9, 8/9, 8/9, 8/9, 8/9, 9/9, 9/9, 9/9, 9/9];
                Amp = containers.Map(key_Amp,value_Amp);
                key_phase =   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
                value_phase = [3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0];
                Phase = containers.Map(key_phase,value_phase);
                for indexes = 1:length(tx_int)
                    tx_symb(indexes) = Amp(tx_int(indexes))*exp(1i*pi*(Phase(tx_int(indexes)) + 1/2*tx_int(indexes)));
                end
                
            case 6 % 16 Orbits
                key_Amp =   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
                value_Amp = [2/17, 2/17, 2/17, 2/17, 3/17, 3/17, 3/17, 3/17, 4/17, 4/17, 4/17, 4/17, 5/17, 5/17, 5/17, 5/17, 6/17, 6/17, 6/17, 6/17, 7/17, 7/17, 7/17, 7/17, 8/17, 8/17, 8/17, 8/17, 9/17, 9/17, 9/17, 9/17, 10/17, 10/17, 10/17, 10/17, 11/17, 11/17, 11/17, 11/17, 12/17, 12/17, 12/17, 12/17, 13/17, 13/17, 13/17, 13/17, 14/17, 14/17, 14/17, 14/17, 15/17, 15/17, 15/17, 15/17, 16/17, 16/17, 16/17, 16/17, 17/17, 17/17, 17/17, 17/17];
                Amp = containers.Map(key_Amp,value_Amp);
                key_phase =   {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63};
                value_phase = [3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0, 3/4, 3/4, 3/4, 3/4, 0, 0, 0, 0];
                Phase = containers.Map(key_phase,value_phase);
                for indexes = 1:length(tx_int)
                    tx_symb(indexes) = Amp(tx_int(indexes))*exp(1i*pi*(Phase(tx_int(indexes)) + 1/2*tx_int(indexes)));
                end
        end
          
    otherwise
        error('You have entered an unsupported Digital Modulation scheme. Please choose among [PSK, PAM, QAM, CQAM].');
end