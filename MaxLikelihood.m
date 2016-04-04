function [rx_symb_ML] = MaxLikelihood(rx_symb,Nbps,modulation)
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 1 - Optimal communication chain over the ideal channel


%% MAXIMUM LIKELIHOOD CRITERION
% By comparing the Received Symbols 'rx_symb' with the Symbols of Reference 'symb_ref' we can apply the ML Criterion and deduce 'rx_symb_ML'
% INPUTS
    % 'rx_symb' - Symbols
    % 'modulation' - Digital Modulation [PSK, PAM, QAM, CQAM]
    % 'Nbps' - Number of bits per symbol [1, 2, 3, 4, 5, 6]
% OUTPUT
    % 'rx_symb_ML' - Symbols after the application of the Maximum Likelihood Criterion

%% INITIALIZATION
Nsymb = 2^Nbps; % Number of Symbols of Reference
switch modulation
    %% ***** PULSE AMPLITUDE MODULATION *****
    case 'PAM'
        % Generating Symbols of Reference
        distance = 2/(2^Nbps-1);
        symb_ref = distance * ((0:2^Nbps-1) - (2^Nbps-1)/2);
        
        % Application of the Maximum Likelihood Criterion
        rx_symb_ML = zeros(1,length(rx_symb)); % Initialize the length of the vector
        for index = 1:length(rx_symb) % At each loop, application of the ML Criterion on a Symbol of the Received Signal
            trans = norm(rx_symb(index) - symb_ref(1)); % Measure the Euclidean Norm with the first Symbol of Reference
            rx_symb_ML(index) = symb_ref(1); % This step is needed in order to apply a comparison from the first inner loop
            
            for index2 = 2:Nsymb % Loop to determine the minimal Euclidean Norm of each Received Symbol
                if norm(rx_symb(index) - symb_ref(index2)) < trans
                    trans = norm(rx_symb(index) - symb_ref(index2));
                    rx_symb_ML(index) = symb_ref(index2);
                end
            end
        end
        rx_symb_ML = rx_symb_ML';
        
        
        
        %% ***** QUADRATURE AMPLITUDE MODULATION *****
    case 'QAM'
        % Generating Symbols of Reference
        Nbps_I = Nbps/2;
        distance_I = 2/(2^Nbps_I-1);
        symb_ref_I = distance_I * ((0:2^Nbps_I-1) - (2^Nbps_I-1)/2);
        Nbps_Q = Nbps/2;
        distance_Q = 2/(2^Nbps_Q-1);
        symb_ref_Q = distance_Q * ((0:2^Nbps_Q-1) - (2^Nbps_Q-1)/2);
        
        t = 1;
        symb_ref = zeros(1,length(symb_ref_I)*length(symb_ref_Q));
        for index = 1:length(symb_ref_I)
            for index2 = 1:length(symb_ref_Q)
                symb_ref(t) = symb_ref_I(index) + 1i*symb_ref_Q(index2);
                t = t + 1;
            end
        end

        % Application of the Maximum Likelihood Criterion
        rx_symb_ML = zeros(1,length(rx_symb));
        for index = 1:length(rx_symb)
            trans = norm(rx_symb(index) - symb_ref(1));
            rx_symb_ML(index) = symb_ref(1);
            
            for index2 = 2:Nsymb % Loop to determine the minimal Euclidean Norm of each Received Symbol
                if norm(rx_symb(index) - symb_ref(index2)) < trans
                    trans = norm(rx_symb(index) - symb_ref(index2));
                    rx_symb_ML(index) = symb_ref(index2);
                end
            end
        end
        rx_symb_ML = rx_symb_ML.';
end
end
