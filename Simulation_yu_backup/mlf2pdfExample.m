%% Script example for mlf2pdf

clear % Clear all variables but not breakpoints
close all % Close all the figures

%% startup.m template
set(0,'defaultFigurePosition', [10 10 600 450]);
set(0,'defaultAxesFontSize', 12);
set(0, 'defaultTextFontSize', 12);
set(0, 'defaultAxesFontName', 'Palatino Linotype');
set(0, 'defaultTextFontName', 'Palatino Linotype');

%% Compute the two functions to plot
x = linspace(-8,8, 31);
yLin = 5*x;
yQuad = x.^2;

figure % Creates a new figure -> avoids erasing the previous one
set(gca, 'FontSize', 13); % Font size for the axes == 13
plot(x, yLin, '-', 'Marker', 's', 'Color', 'black',...
    'MarkerSize', 6, 'LineWidth', 1,...
    'DisplayName', '$y(x) = 5x$');
hold on; % Allows to plot several curves on the same figure
plot(x, yQuad, '-', 'Marker', 'o', 'Color', 'blue',...
    'MarkerSize', 6, 'LineWidth', 2,...
    'DisplayName', '$y(x) = x^2$');
hold off; % Disables hold on
grid on; % Shows the grid
xlabel('Value of $x$');
ylabel('Value of $y(x)$');
title('A plot of two functions');

%% Other possibly interesting commands:
% axis tight % Make the axis span thighly the extension of the curves
% xlim([xMin, xMax]); % Show only the part of the x axis going from xMin to
                    % xMax
% ylim([yMin, yMax]); % Show only the part of the y axis going from yMin to
                    % yMax
% set(gca,'YDir','reverse'); % Reverses the orientation of the y axis
% set(gca,'YDir','normal'); % Makes the orientation of the y axis standard

%% Show the legend
% Show the legend using the property 'DisplayName' of each curve that has
% been plotted
hLeg = legend('show');

% Modify the font size for the legend. Set LaTeX as interpreter for the
% display and the computation of the legend box boundaries. Position works
% like this: [(bottom left corner horizontal position), (bottom left corner
% vertical position), (width), (height)]. If the width and/or height are
% too small for proper display, then Matlab adjusts the parameters => 
% interesting to set 0.01 to have the tighest possible legend box.
set(hLeg, 'fontsize', 12, 'interpreter', 'latex',...
    'Position', [0.76,0.2,0.01,0.01]);

%% Save to PDF
% Now, we can use mlf2pdf to create the figure. The figure is Figure 1 as
% only one figure has been plotted and the possibly preexisting figures
% have been closed by executing 'close all'.
% mlf2pdf(1, 'A_Figure');
