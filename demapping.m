function [rx_bin] = demapping(rx_symb,Nbps,modulation)
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 1 - Optimal communication chain over the ideal channel


%% DEMAPPING
% INPUTS
    % 'rx_symb' - Symbols
    % 'modulation' - Digital Modulation [PSK, PAM, QAM, CQAM]
    % 'Nbps' - Number of bits per symbol [1, 2, 3, 4, 5, 6]
% OUTPUT
    % 'rx_bit' - Bitstream

%% INITIALIZATION
Nsymb = size(rx_symb,1); % Number of symbols

switch modulation
    %% ***** PULSE AMPLITUDE MODULATION *****
    case 'PAM'
        % Symbol to integer
        distance = 2/(2^Nbps-1);
        rx_int = rx_symb/distance + (2^Nbps-1)/2;
        % Integer to binary
        rx_int = real(round(rx_int));
        rx_int = max(rx_int, 0); % Discard the potential values interpreted as smaller than the smallest possible integer (0)
        rx_int(rx_int > 2^Nbps - 1) = 2^Nbps - 1; % Discard the potential values interpreted as bigger than the biggest possible integer (2^Nbps-1)
        rx_bin  = fliplr(de2bi(rx_int));
        % Binary to gray
        for ii = 2:Nbps,
            rx_bin(:,ii) = xor( rx_bin(:,ii-1) , rx_bin(:,ii) );
        end
        rx_bin = reshape(rx_bin',Nsymb*Nbps,1);

        
        
    %% ***** QUADRATURE AMPLITUDE MODULATION *****    
    case  'QAM'
       switch Nbps
            case 1 % Go back to the simple PAM case
                distance = 2/(2^Nbps-1);
                rx_int = rx_symb/distance + (2^Nbps-1)/2;
                rx_int = real(round(rx_int));
                rx_int = max(rx_int, 0);
                rx_int(rx_int > 2^Nbps - 1) = 2^Nbps - 1;
                rx_bin  = fliplr(de2bi(rx_int));
                for ii = 2:Nbps,
                    rx_bin(:,ii) = xor( rx_bin(:,ii-1) , rx_bin(:,ii) );
                end
                rx_bin = reshape(rx_bin',Nsymb*Nbps,1);

            otherwise % RECTANGULAR-QAM
        
            Nbps_I = ceil(Nbps/2); % Number of symbols carried in the 'Inphase Carrier'
            Nbps_Q = Nbps - Nbps_I; % Number of symbols carried in the 'Quadrature Carrier'

            %% *** REAL PART ***
            rx_symb_I = real(rx_symb);
            % Symbol to integer
            distance = 2/(2^Nbps_I-1); 
            rx_int_I = rx_symb_I/distance + (2^Nbps_I-1)/2;
            % Integer to binary
            rx_int_I = round(rx_int_I);
            rx_int_I = max(rx_int_I, 0); % Discard the potential values interpreted as smaller than the smallest possible integer (0)
            rx_int_I(rx_int_I > 2^Nbps_I - 1) = 2^Nbps_I - 1; % Discard the potential values interpreted as bigger than the biggest possible integer (2^Nbps-1)
            rx_bin_I  = fliplr(de2bi(rx_int_I));
            % Binary to gray
            for ii = 2:Nbps_I,
                rx_bin_I(:,ii) = xor( rx_bin_I(:,ii-1) , rx_bin_I(:,ii) );
            end

            %% *** IMAGINARY PART ***
            rx_symb_Q = imag(rx_symb);
            % Symbol to integer
            rx_int_Q = rx_symb_Q/distance + (2^Nbps_Q-1)/2;
            % Integer to binary
            rx_int_Q = round(rx_int_Q);
            rx_int_Q = max(rx_int_Q, 0); % Discard the potential values interpreted as smaller than the smallest possible integer (0)
            rx_int_Q(rx_int_Q > 2^Nbps_Q - 1) = 2^Nbps_Q - 1; % Discard the potential values interpreted as bigger than the biggest possible integer (2^Nbps-1)
            rx_bin_Q  = fliplr(de2bi(rx_int_Q));
            % Binary to gray
            for ii = 2:Nbps_Q,
                rx_bin_Q(:,ii) = xor( rx_bin_Q(:,ii-1) , rx_bin_Q(:,ii) );
            end     

            %% *** BIT CONCATENATION ***
            rx_bin = reshape([rx_bin_I,rx_bin_Q]',Nsymb*Nbps,1);
       end
        
        
    %% ***** PHASE SHIFT KEYING *****
    case 'PSK'
        % Symbol to integer
        rx_int = (pi*5/4 + angle(rx_symb)) * (2^Nbps)/(pi*2);
        % Integer to binary
        rx_int = real(round(rx_int));
        rx_int = max(rx_int, 0); % Discard the potential values interpreted as smaller than the smallest possible integer (0)
        rx_int(rx_int > 2^Nbps - 1) = 2^Nbps - 1; % Discard the potential values interpreted as bigger than the biggest possible integer (2^Nbps-1)
        rx_bin  = fliplr(de2bi(rx_int));
        % Binary to gray
        for ii = 2:Nbps,
            rx_bin(:,ii) = xor( rx_bin(:,ii-1) , rx_bin(:,ii) );
        end
        rx_bin = reshape(rx_bin',Nsymb*Nbps,1);
end