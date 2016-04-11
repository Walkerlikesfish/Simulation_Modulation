function FoM()
% Modulation and Coding Project
% TEAM: MOY - Mroueh Michael, Asfour A. Omar, Liu Yu
% April 2016
% Part 1 - Optimal communication chain over the ideal channel


%% This function is used in order to plot different Graphs and FoMs
clear; clc; close all;


%% DECOMMENTS HERE WHAT YOU WANT TO PLOT 'print -dmeta' command to add to .docx
% rx_Constellation_Symbols('PAM', 3, 10); % Plot the Transmitted Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('PSK', 1, 3); % Plot the Transmitted Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('Cross-QAM', 3, 10); % Plot the Transmitted Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('Optimal-8QAM', 3, 10); % Plot the Transmitted Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('QAM', 4, 10); % Plot the Transmitted Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('PAM', 4, 10); % Plot the Transmitted Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('PSK', 5, 10); % Plot the Received Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('QAM', 4, 10); % Plot the Received Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('Cross-QAM', 7, 10); % Plot the Received Constellation Symbols (1 Graph)
% rx_Constellation_Symbols('Optimal-8QAM', 3, 10); % Plot the Received Constellation Symbols (1 Graph)

% BER_EbN0_NGraphs(); % Plot BER in function of Eb/N0 for different Modulation Schemes on [2, 4, 8, 16, 32, 64,...] Symbols (N Graphs)
% BER_EbN0_1Graph('QAM'); % Plot BER in function of Eb/N0 for the Modulation Scheme given in argument ('QAM' or 'PAM' or 'Cross-QAM' or 'PSK')' on [2, 4, 8, 16, 32, 64] Symbols (1 Graph)
% SER_EbN0_NGraphs(); % Plot SER in function of Eb/N0 for different Modulation Schemes on [2, 4, 8, 16, 32, 64,...] Symbols (N Graphs)
% Hamming_Distance(); % Plot Ecoding ('BER/SER*Nbps') in function of Eb/N0 for different Modulation Schemes on [2, 4, 8, 16, 32, 64,...] Symbols (N Graphs) - Corresponds to a 'Pseudo Apparent Hamming Distance' 

%tx_Constellation_Symbols('QAM', 3, -1)

% RRCF_IR_Spectrum(0.25, 33, 1/(5e6), 4); % Plot the Spectrum and Impulse Response of the RRCF for different Roll-off Factors and for a given Number of Taps

 

%% ********************
%% ** PLOT FUNCTIONS **
%% ********************

%% Plot the Transmitted Constellation Symbols (1 Graph)
function tx_Constellation_Symbols(mod, Nbps, EbN0)
[~, ~, tx_symb, ~] = DVBS2CommunicationChain(mod, Nbps, EbN0);
figure; scatter(real(tx_symb),imag(tx_symb), 'filled');
ax=gca; ax.XAxisLocation = 'bottom'; ax.YAxisLocation = 'left'; % Set the axis to the origin
ax.XLim = [min(real(tx_symb)) - 0.25, max(real(tx_symb)) + 0.25]; ax.YLim = [min(real(tx_symb)) - 0.25, max(real(tx_symb)) + 0.25]; ax.XLimMode = 'manual'; axis square; % Set the axis limits
title(['Constellation Symbols of ', num2str(2^Nbps), mod]); xlabel('In-Phase'); ylabel('Quadrature'); hold on;


function rx_Constellation_Symbols(mod, Nbps, EbN0)
%% Plot the Received Constellation Symbols (1 Graph)
[~, ~, tx_symb, rx_symb] = DVBS2CommunicationChain(mod, Nbps, EbN0);
figure; scatter(real(rx_symb),imag(rx_symb), 'magenta', 'filled'); hold on; scatter(real(tx_symb),imag(tx_symb), 20, 'blue', '+');
ax=gca; ax.XAxisLocation = 'bottom'; ax.YAxisLocation = 'left';
ax.XLim = [min(real(tx_symb)) - 0.25, max(real(tx_symb)) + 0.25]; ax.YLim = [min(real(tx_symb)) - 0.25, max(real(tx_symb)) + 0.25]; ax.XLimMode = 'manual'; axis square;
title(['Constellation Symbols of ', num2str(2^Nbps), mod, ' with E_{b}/N_{0} = ', num2str(EbN0)]); xlabel('In-Phase'); ylabel('Quadrature');


%% Plot BER in function of Eb/N0 for different Modulation Schemes on [2, 4, 8, 16, 32, 64,...] Symbols (N Graphs)
function BER_EbN0_NGraphs()
min_Eb_N0 = 0; % Minimal value given to the Eb/N0 Ratio
step_Eb_N0 = 0.40;
max_Eb_N0_12 = 4; % Maximal value given to the Eb/N0 Ratio for 2 and 4 Symbols (Nbps = [1, 2])
max_Eb_N0 = max_Eb_N0_12; % Value set to 'max_Eb_N0_12' for 2 and 4 Symbols. Value increased to 'max_Eb_N0_12*[2, 3, 4, 5,...]' for [8, 16, 32, 64,...] symbols
% This is done in order to compute relevant points
BER = zeros(5,(round((max_Eb_N0*5-min_Eb_N0)/step_Eb_N0) + 1));

for Nbps = 1:6
    index = 1;
    max_Eb_N0 = max_Eb_N0_12 * (Nbps);
    for EbN0 = min_Eb_N0:step_Eb_N0:max_Eb_N0
        [BER(1,index), ~, ~, ~] = DVBS2CommunicationChain('PAM', Nbps, EbN0);
        [BER(2,index), ~, ~, ~] = DVBS2CommunicationChain('PSK', Nbps, EbN0);
        [BER(3,index), ~, ~, ~] = DVBS2CommunicationChain('Cross-QAM', Nbps, EbN0);
        [BER(4,index), ~, ~, ~] = DVBS2CommunicationChain('QAM', Nbps, EbN0);
        if Nbps == 3
            [BER(5,index), ~, ~, ~] = DVBS2CommunicationChain('Optimal-8QAM', Nbps, EbN0);
        end
        index = index + 1;
    end
    
    figure;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(1, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(2, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(3, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(4, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
    title({('Figure of Merit of the BER in function of the E_{b}/N_{0} ratio'), ['Comparison of different Modulation Schemes for ' num2str(2^Nbps), ' symbols']}); xlabel('E_{b}/N_{0} [dB]'); ylabel('BER');
    legend([num2str(2^Nbps),'PAM'], [num2str(2^Nbps),'PSK'], [num2str(2^Nbps),'Cross-QAM'], [num2str(2^Nbps),'QAM'], 'Location','southoutside'); legend('boxoff');
    if Nbps == 3
        semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(5, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
        legend([num2str(2^Nbps),'PAM'], [num2str(2^Nbps),'PSK'], [num2str(2^Nbps),'Cross-QAM'], [num2str(2^Nbps),'QAM'], 'Optimal-8QAM', 'Location','southoutside'); legend('boxoff');
    end
end


%% Plot BER in function of Eb/N0 for a Modulation Scheme ('QAM' or 'PAM' or 'Cross-QAM' or 'PSK')' on [2, 4, 8, 16, 32, 64] Symbols (1 Graph)
function BER_EbN0_1Graph(mod)
min_Eb_N0 = 0; % Minimal value given to the Eb/N0 Ratio
step_Eb_N0 = 0.4;
max_Eb_N0 = 10; % Maximal value given to the Eb/N0 Ratio
BER = zeros(1,(round((max_Eb_N0-min_Eb_N0)/step_Eb_N0) + 1));

figure;
for Nbps = 1:6
    index = 1;
    for EbN0 = min_Eb_N0:step_Eb_N0:max_Eb_N0
        [BER(1,index), ~, ~, ~] = DVBS2CommunicationChain(mod, Nbps, EbN0);
        index = index + 1;
    end
    
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(1, :), '-o'); hold on;
end
title({('Figure of Merit of the BER in function of the E_{b}/N_{0} ratio'), ['Comparison of the ', mod, ' for [2, 4, 8, 16, 32, 64] symbols']}); xlabel('E_{b}/N_{0} [dB]'); ylabel('BER');
legend(['2', mod], ['4', mod], ['8', mod], ['16', mod], ['32', mod], ['64', mod], 'Location','southoutside'); legend('boxoff');


%% Plot SER in function of Eb/N0 for different Modulation Schemes on [2, 4, 8, 16, 32, 64,...] Symbols (N Graphs)
function SER_EbN0_NGraphs()
min_Eb_N0 = 0; % Minimal value given to the Eb/N0 Ratio
step_Eb_N0 = 0.40;
max_Eb_N0_12 = 4;
max_Eb_N0 = max_Eb_N0_12;
SER = zeros(5,(round((max_Eb_N0*5-min_Eb_N0)/step_Eb_N0) + 1));

for Nbps = 1:6
    index = 1;
    max_Eb_N0 = max_Eb_N0_12 * (Nbps);
    for EbN0 = min_Eb_N0:step_Eb_N0:max_Eb_N0
        [~, SER(1,index), ~, ~] = DVBS2CommunicationChain('PAM', Nbps, EbN0);
        [~, SER(2,index), ~, ~] = DVBS2CommunicationChain('PSK', Nbps, EbN0);
        [~, SER(3,index), ~, ~] = DVBS2CommunicationChain('Cross-QAM', Nbps, EbN0);
        [~, SER(4,index), ~, ~] = DVBS2CommunicationChain('QAM', Nbps, EbN0);
        if Nbps == 3
            [~, SER(5,index), ~, ~] = DVBS2CommunicationChain('Optimal-8QAM', Nbps, EbN0);
        end
        index = index + 1;
    end
    
    figure;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, SER(1, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, SER(2, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, SER(3, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, SER(4, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
    title({('Figure of Merit of the SER in function of the E_{b}/N_{0} ratio'), ['Comparison of different Modulation Schemes for ' num2str(2^Nbps), ' symbols']}); xlabel('E_{b}/N_{0} [dB]'); ylabel('SER');
    legend([num2str(2^Nbps),'PAM'], [num2str(2^Nbps),'PSK'], [num2str(2^Nbps),'Cross-QAM'], [num2str(2^Nbps),'QAM'], 'Location','southoutside'); legend('boxoff');
    if Nbps == 3
        semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, SER(5, 1:round(max_Eb_N0/step_Eb_N0)+1), '-o'); hold on;
        legend([num2str(2^Nbps),'PAM'], [num2str(2^Nbps),'PSK'], [num2str(2^Nbps),'Cross-QAM'], [num2str(2^Nbps),'QAM'], 'Optimal-8QAM', 'Location','southoutside'); legend('boxoff');
    end
end


%% Plot Ecoding ('BER/SER*Nbps') in function of Eb/N0 for different Modulation Schemes on [2, 4, 8, 16, 32, 64,...] Symbols (N Graphs) - Corresponds to a 'Pseudo Apparent Hamming Distance' 
function Hamming_Distance()
min_Eb_N0 = 0;
step_Eb_N0 = 0.40;
max_Eb_N0_12 = 4;
max_Eb_N0 = max_Eb_N0_12;
BER = zeros(5,(round((max_Eb_N0*5-min_Eb_N0)/step_Eb_N0) + 1));
SER = zeros(5,(round((max_Eb_N0*5-min_Eb_N0)/step_Eb_N0) + 1));

for Nbps = 1:6
    index = 1;
    max_Eb_N0 = max_Eb_N0_12 * (Nbps);
    for EbN0 = min_Eb_N0:step_Eb_N0:max_Eb_N0
        [BER(1,index), SER(1,index), ~, ~] = DVBS2CommunicationChain('PAM', Nbps, EbN0);
        [BER(2,index), SER(2,index), ~, ~] = DVBS2CommunicationChain('PSK', Nbps, EbN0);
        [BER(3,index), SER(3,index), ~, ~] = DVBS2CommunicationChain('Cross-QAM', Nbps, EbN0);
        [BER(4,index), SER(4,index), ~, ~] = DVBS2CommunicationChain('QAM', Nbps, EbN0);
        if Nbps == 3
            [BER(5,index), SER(5,index), ~, ~] = DVBS2CommunicationChain('Optimal-8QAM', Nbps, EbN0);
        end
        index = index + 1;
    end
    
    figure;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(1, 1:round(max_Eb_N0/step_Eb_N0)+1)./SER(1, 1:round(max_Eb_N0/step_Eb_N0)+1)*Nbps, '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(2, 1:round(max_Eb_N0/step_Eb_N0)+1)./SER(2, 1:round(max_Eb_N0/step_Eb_N0)+1)*Nbps, '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(3, 1:round(max_Eb_N0/step_Eb_N0)+1)./SER(3, 1:round(max_Eb_N0/step_Eb_N0)+1)*Nbps, '-o'); hold on;
    semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(4, 1:round(max_Eb_N0/step_Eb_N0)+1)./SER(4, 1:round(max_Eb_N0/step_Eb_N0)+1)*Nbps, '-o'); hold on;
    title({('Apparent Hamming Distance \eta in function of the E_{b}/N_{0} ratio '), ['Comparison of different Modulation Schemes for ' num2str(2^Nbps), ' symbols'], ' \eta = BER/SER*Nbps'}); xlabel('E_{b}/N_{0} [dB]'); ylabel('\eta');
    legend([num2str(2^Nbps),'PAM'], [num2str(2^Nbps),'PSK'], [num2str(2^Nbps),'Cross-QAM'], [num2str(2^Nbps),'QAM'], 'Location','southoutside'); legend('boxoff');
    if Nbps == 3
        semilogy(min_Eb_N0:step_Eb_N0:max_Eb_N0, BER(5, 1:round(max_Eb_N0/step_Eb_N0)+1)./SER(5, 1:round(max_Eb_N0/step_Eb_N0)+1)*Nbps, '-o'); hold on;
        legend([num2str(2^Nbps),'PAM'], [num2str(2^Nbps),'PSK'], [num2str(2^Nbps),'Cross-QAM'], [num2str(2^Nbps),'QAM'], 'Optimal-8QAM', 'Location','southoutside'); legend('boxoff');
    end
end


%% Plot the Spectrum and Impulse Response of the RRCF for different Roll-off Factors and Number of Taps
function RRCF_IR_Spectrum(Beta_step, Ntaps, Tsymb, M)
fig1 = figure;
fig2 = figure;
for Beta = 0:Beta_step:1
    [IR_RRC, H_RRCF, t_axis, f_axis] = RRCFDesign(Beta, Ntaps, Tsymb, M/Tsymb);
    figure(fig1); plot(t_axis, IR_RRC, '-o'); hold on;
    figure(fig2); plot(f_axis, H_RRCF, '-o'); hold on;
end
legendCell = cellstr(num2str((0:Beta_step:1)', 'Beta = %g'));
figure(fig1); title({'Impulse Response of the Root-Raised-Cosine Filter', ['Comparison for different values of the Roll-off Factor with N_{taps} = ' num2str(Ntaps) ', T_{symb} = ' num2str(Tsymb) ' and M = ' num2str(M)]}); xlabel('Time [µs]'); ylabel('g(t)'); legend(legendCell, 'Location','southoutside'); legend('boxoff');
figure(fig2); title({'Spectrum of the Root-Raised-Cosine Filter', ['Comparison for different values of the Roll-off Factor with N_{taps} = ' num2str(Ntaps) ', T_{symb} = ' num2str(Tsymb) ' and M = ' num2str(M)]}); xlabel('Frequency [MHz]'); ylabel('G(f)'); legend(legendCell, 'Location','southoutside'); legend('boxoff');
