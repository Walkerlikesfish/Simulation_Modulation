function [rx_symb_ML] = maxLikelihood(rx_symb,Nbps,modulation)
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 1 - Optimal communication chain over the ideal channel


%% MAXIMUM LIKELIHOOD CRITERION
% By comparing the Received Symbols 'rx_symb' with the Symbols of Reference 'symb_ref' we can apply the ML Criterion and deduce 'rx_symb_ML'
% INPUTS
%       'rx_symb' - Symbols
%       'modulation' - Digital Modulation [PSK, PAM, QAM, Cross-QAM, Optimal-8QAM]
%       'Nbps' - Number of bits per symbol
% OUTPUT
%       'rx_symb_ML' - Symbols after the application of the Maximum Likelihood Criterion


%% ALGORITHM
% 'symb_ref' is a [2^Nbps x 2] matrix used as a Containers.map
% The 1st column includes all the Symbols of Reference in their Complex form (-0.7071 + i*0.7071, +0.7071 - i*0.7071,...)
% The 2nd column includes the associated Decimal Value (0, 1, 2, 3, ...)
% Step 1: Generating Symbols of Reference in function of (Nbps, modulation)
% Step 2: Application of the Maximum Likelihood Criterion -> Compare each element of 'rx_symb' with the Symbols of Reference 'symb_ref'
% Step 3: Return the Symbols after the application of the Maximum Likelihood Criterion 'rx_symb_ML'

switch modulation
        %% ***** PULSE AMPLITUDE MODULATION *****
    case 'PAM'
        % Generating Symbols of Reference
        symb_ref = zeros(2^Nbps,2);
        distance = 2/(2^Nbps-1);
        for index = 1:2^Nbps
            symb_ref(index,:) = [distance*(index-1 - (2^Nbps-1)/2), index-1];
        end
        % Application of the Maximum Likelihood Criterion
        rx_symb_ML = appML(rx_symb, Nbps, symb_ref);
        
        
        %% ***** QUADRATURE AMPLITUDE MODULATION *****
    case 'QAM'
        % Generating Symbols of Reference
        Nbps_I = ceil(Nbps/2); Nbps_Q = Nbps - Nbps_I;
        symb_ref = zeros(2^Nbps,2);
        distance = 2/(2^Nbps_I-1);
        for index_I = 1:2^Nbps_I
            for index_Q = 1:2^Nbps_Q
                symb_ref(index_I + (2^Nbps_I)*(index_Q-1), :) = [distance*((index_I-1 - (2^Nbps_I-1)/2) - 1i* (index_Q-1 - (2^Nbps_Q-1)/2)), index_I-1 + (2^Nbps_I)*(index_Q-1)];
            end
        end
        % Application of the Maximum Likelihood Criterion
        rx_symb_ML = appML(rx_symb, Nbps, symb_ref);
        
        %% ***** PHASE SHIFT KEYING *****
    case 'PSK'
        symb_ref = zeros(2^Nbps,2);
        % Generating Symbols of Reference
        for index = 1:2^Nbps
            symb_ref(index,:) = [exp(1i*pi*(3/4 + 2/(2^Nbps)*(index-1))), index-1];
        end
        % Application of the Maximum Likelihood Criterion
        rx_symb_ML = appML(rx_symb, Nbps, symb_ref);
        
        %% ***** CROSS-QUADRATURE AMPLITUDE MODULATION *****
    case 'Cross-QAM'
        % Generating Symbols of Reference
        symb_ref = zeros(2^Nbps,2);
        switch Nbps
            case 1 % Only 2 Symbols on the 1st Orbit => Elementary PSK Case
                symb_ref = [exp(1i*pi*(3/4 + 2/(2^Nbps)*0)), 0; exp(1i*pi*(3/4 + 2/(2^Nbps)*1)), 1];
            otherwise
                for index = 1:2^Nbps
                    radius = (2 + floor((index - 1)/4)) / (2^Nbps/4 + 1);
                    phase = (mod(floor((index-1)/4), 2) == 0)*3/4;
                    symb_ref(index,:) = [radius*exp(1i*pi*(phase + 1/2*mod((index - 1),4))), index - 1];
                end
        end
        % Application of the Maximum Likelihood Criterion
        rx_symb_ML = appML(rx_symb, Nbps, symb_ref);
        
        
        %% ***** OPTIMAL-8QAM *****
    case 'Optimal-8QAM'
        % Generating Symbols of Reference
        symb_ref = zeros(2^Nbps,2);
        radius = [sqrt(2)/(1 + sqrt(3)), 1];
        for index = 1:2^Nbps
            phase = (floor((index-1)/4) == 0)*3/4;
            symb_ref(index,:) = [radius(floor((index-1)/4) + 1)*exp(1i*pi*(phase + 1/2*mod((index - 1),4))), index - 1];
        end
        % Application of the Maximum Likelihood Criterion
        rx_symb_ML = appML(rx_symb, Nbps, symb_ref);
end


%% ***** APPLICATION OF THE MAXIMUM LIKELIHOOD CRITERION *****
function [rx_symb_ML] = appML(rx_symb, Nbps, symb_ref)
rx_symb_ML = zeros(length(rx_symb), 1); % Initialize the size of 'rx_symb_ML' in order to save some time
CSM = sqrt(abs((rx_symb(:) - symb_ref(1, 1)) .^ 2)); % CSM - Current Smallest Norm: Transient save allowing to 'compare before assign'

for index = 1:2^Nbps % Loop to determine the minimal Euclidean Norm of each Received Symbol
    BOOL = CSM > sqrt(abs((rx_symb(:) - symb_ref(index, 1)) .^ 2)); % Boolean Vector taking a true value ('1') when the corresponding CSM row appears to not be the smallest possible Euclidean Norm
    CSM(BOOL) = sqrt(abs((rx_symb(BOOL) - symb_ref(index, 1)) .^ 2)); % Update CSM accordingly
    rx_symb_ML(BOOL) = symb_ref(index, 2);
end