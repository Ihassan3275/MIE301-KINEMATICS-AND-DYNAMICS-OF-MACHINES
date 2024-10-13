%%% MIE301 Lab 1 - Case 1: Link 3 hangs free under gravity
%% This file contains code from MIE301 instructors and Ibrahim Hassan
%% This file includes all MATLAB code used for parts (a), (b), and (c)
%% To generate expected results, you will have to comment/uncomment the lines of code for each part.
%%
close all; % closes all figures
clear all; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace

% Plot Parameters: these will be used to set the axis limits on the figure
xmin = -15;  % leftmost window edge
xmax = 15;   % rightmost window edge
ymin = -20;  % bottom window edge
ymax = 15;   % top window edge

%% Define the link and motion Parameters (part (a))
%stepsize = 6 *pi/180;                        % step size
%max_rotation_theta2 = 320 *pi/180;            % theta2 from 0 to 320 degrees
%theta2 = 0 : stepsize : max_rotation_theta2;  % rotation steps of link 2

%% Define link and motion Parameters (part (b))
%theta2 = [0, 80, 160, 240, 320] * pi / 180;  % specific configurations for stroboscopic effect (in radians)
%theta_labels = {'0°', '80°', '160°', '240°', '320°'};  % Labels for the configurations

%% Part (c) [N = 6, 10, 20]:
%% Define theta2 for N=6 (0°, 64°, 128°, 192°, 256°, 320° in radians)
N = 6;
theta2 = linspace(0, 320 * pi / 180, N);  % N equally spaced angles
theta_labels = {'0°', '64°', '128°', '192°', '256°', '320°'};

%% Define theta2 for N=10 (equally spaced in radians from 0° to 320°)
% N = 10;
% theta2 = linspace(0, 320 * pi / 180, N);  % N equally spaced angles
% theta_labels = arrayfun(@(x) [num2str(x), '°'], linspace(0, 320, N), 'UniformOutput', false);  % Generates labels for angles

%% Define theta2 for N=20 (equally spaced in radians from 0° to 320°)
% N = 20;
% theta2 = linspace(0, 320 * pi / 180, N);  % N equally spaced angles
% theta_labels = arrayfun(@(x) [num2str(x), '°'], linspace(0, 320, N), 'UniformOutput', false);  % Generates labels for angles

%% Case 1 - Link lengths
length2 = 8;    % link 2 length in case 1
length3 = 4;    % link 3 length in case 1

%% Plot Case 1 Figure
figure(1);                         % create new figure
set(1, 'WindowStyle', 'Docked');   % dock the figure

% Initialize arrays for path storage
Bx_path = zeros(1, length(theta2));
By_path = zeros(1, length(theta2));
Cx_path = zeros(1, length(theta2));
Cy_path = zeros(1, length(theta2));

%% Calculate mechanism motion and plot it
for i = 1:length(theta2)                     % step through the 5 specified configurations
    % Point A (pivot)
    Ax(i) = 0;                                % pivot point of link 2 position
    Ay(i) = 0;                                % pivot point of link 2 position
  
    % Point B (end of link 2)
    Bx(i) = length2 * cos(theta2(i));         % point B position
    By(i) = length2 * sin(theta2(i)); 
    
    % Point C (end of link 3)
    Cx(i) = Bx(i);                            % Corrected index usage for Cx
    Cy(i) = By(i) - length3;                  % Corrected index usage for Cy

    % Store positions of B and C for path plotting later
    Bx_path(i) = Bx(i);
    By_path(i) = By(i);
    Cx_path(i) = Cx(i);
    Cy_path(i) = Cy(i);
    
    % Plot links for the current configuration
    plot([Ax(i) Bx(i)], [Ay(i) By(i)], 'r', 'LineWidth', 3);  % Link 2
    hold on;
    plot([Bx(i) Cx(i)], [By(i) Cy(i)], 'b', 'LineWidth', 2);  % Link 3
    
    % Plot points
    plot(Bx(i), By(i), 'bo', 'MarkerFaceColor', 'w');         % Point B
    plot(Cx(i), Cy(i), 'go', 'MarkerFaceColor', 'w');         % Point C
    
    %% Label points
    text(Bx(i) + 0.4, By(i), 'B', 'color', 'b');
    text(Cx(i) + 0.4, Cy(i), 'C', 'color', 'g');
    
    %% Label configurations with theta values
    text(Bx(i) + 1, By(i) + 1, ['\theta_2 = ', theta_labels{i}], 'FontSize', 10, 'color', 'k');
    
    %% Draw mass at point C as a box
    rectangle('Position', [Cx(i)-0.6, Cy(i)-0.6, 1.0, 1.0], 'FaceColor', 'g');
    text(Cx(i) - 0.4, Cy(i), 'm', 'FontSize', 5, 'Color', [0 0 0]); % adding the 'm' (for mass) label to the rectangle

    %% Set axis and labels
    axis([xmin xmax ymin ymax]);    % set figure axis limits
    axis equal;                     % ensure equal scaling for x and y
    xlabel('x (cm)', 'fontsize', 15);  % x-axis label
    ylabel('y (cm)', 'fontsize', 15);  % y-axis label
    grid on;                         % add a grid to the figure

    %% Part (a) - title:
    % title('Lab 1 - Case 1');

    %% Part (b) - title:
    % title('Lab 1 - Case 1 with 5 Configurations');

    %% Part (c1) - title:
    title('Lab 1 - Case 1 with N=6');

    %% Part (c2) - title:
    % title('Lab 1 - Case 1 with N=10');

    %% Part (c3) - title:
    % title('Lab 1 - Case 1 with N=20');
end

% Plot the path of points B and C
plot(Bx_path, By_path, '--g', 'linewidth', 1);  % Path of point B
plot(Cx_path, Cy_path, '--b', 'linewidth', 1);  % Path of point C

% Final configuration: hold the last plotted configuration for display
plot([Ax(end) Bx_path(end)], [Ay(end) By_path(end)], 'r', 'LineWidth', 1);  % Final position of Link 2
plot([Bx_path(end) Cx_path(end)], [By_path(end) Cy_path(end)], 'b', 'LineWidth', 1);  % Final position of Link 3
plot(Bx_path(end), By_path(end), 'bo', 'MarkerFaceColor', 'w');   % Final Point B
plot(Cx_path(end), Cy_path(end), 'go', 'MarkerFaceColor', 'w');   % Final Point C

% Show the figure
hold off;