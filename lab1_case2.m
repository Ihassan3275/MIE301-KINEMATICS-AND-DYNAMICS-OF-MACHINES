%%% MIE301 Lab 1 - Case 2: Rigid connection between links 2 and 3
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
% stepsize = 6 *pi/180;                         % step size
% max_rotation_theta2 = 320 *pi/180;            % theta2 from 0 to 320 degrees
% theta2 = 0 : stepsize : max_rotation_theta2;  % rotation steps of link 2

%% Define link and motion Parameters (part (b))
% stepsize = 80 * pi / 180;                              % step size for stroboscopic effect (80 degrees in radians)
% theta2 = [0, 80, 160, 240, 320] * pi / 180;            % specific configurations to plot (in radians)
% theta_labels = {'0°', '80°', '160°', '240°', '320°'};  % Labels for the configurations

%% Part (c) [N = 6, 10, 20]:
%% Define theta2 for N=6 (0°, 64°, 128°, 192°, 256°, 320° in radians)
% N = 6;
% theta2 = linspace(0, 320 * pi / 180, N);  % N equally spaced angles
% theta_labels = {'0°', '64°', '128°', '192°', '256°', '320°'};

%% Define theta2 for N=10 (equally spaced in radians from 0° to 320°)
% N = 10;
% theta2 = linspace(0, 320 * pi / 180, N);  % N equally spaced angles
% theta_labels = arrayfun(@(x) [num2str(x), '°'], linspace(0, 320, N), 'UniformOutput', false); % Generates labels for angles

%% Define theta2 for N=20 (equally spaced in radians from 0° to 320°)
N = 20;
theta2 = linspace(0, 320 * pi / 180, N);  % N equally spaced angles
theta_labels = arrayfun(@(x) [num2str(x), '°'], linspace(0, 320, N), 'UniformOutput', false); % Generates labels for angles

%% Case 2 - Link lengths
length2 = 4;    % link 2 length in case 2
length3 = 8;    % link 3 length in case 2

%% Constant angle phi
phi = 69 * pi / 180;    % angle in radians (69 degrees, last 2 digits of student number)

% Initialize arrays to store the path of points B and C
Bx_path = zeros(1, length(theta2));
By_path = zeros(1, length(theta2));
Cx_path = zeros(1, length(theta2));
Cy_path = zeros(1, length(theta2));

% Plot Case 2 - Rigid connection
figure(1);
set(1,'WindowStyle','Docked')

for i = 1:length(theta2)
    % Point A (pivot)
    Ax = 0; 
    Ay = 0;
    
    % Point B (end of link 2)
    Bx = length2 * cos(theta2(i));   
    By = length2 * sin(theta2(i));   
    
    % Point C (end of link 3, rigid connection)
    Cx = Bx + length3 * cos(theta2(i) + phi);  
    Cy = By + length3 * sin(theta2(i) + phi);
    
    % Store the current positions of points B and C
    Bx_path(i) = Bx;
    By_path(i) = By;
    Cx_path(i) = Cx;
    Cy_path(i) = Cy;
    
    % Plot links for the current configuration
    plot([Ax Bx], [Ay By], 'r', 'LineWidth', 3);  % Link 2
    hold on;
    plot([Bx Cx], [By Cy], 'b', 'LineWidth', 2);  % Link 3
    
    % Plot points
    plot(Bx, By, 'bo', 'MarkerFaceColor', 'w');   % Point B
    plot(Cx, Cy, 'go', 'MarkerFaceColor', 'w');   % Point C
    
     %% Label points
    text(Bx + 0.4, By, 'B', 'color', 'b');
    text(Cx + 0.4, Cy, 'C', 'color', 'g');
    
    %% Label configurations with theta values
    text(Bx + 1, By + 1, ['\theta_2 = ', theta_labels{i}], 'FontSize', 10, 'color', 'k');

    %% Set axis and labels
    axis([xmin xmax ymin ymax]);
    axis equal;
    grid on;
    xlabel('x (cm)');
    ylabel('y (cm)');
    %% Part (a) - title:
    % title('Lab 1 - Case 2');

    %% Part (b) - title:
    % title('Lab 1 - Case 2 with 5 Configurations');

    %% Part (c1) - title:
    % title('Lab 1 - Case 2 with N=6');

    %% Part (c2) - title:
    % title('Lab 1 - Case 2 with N=10');

    %% Part (c3) - title:
    title('Lab 1 - Case 2 with N=20');
end

% Plot the path of points B and C
plot(Bx_path, By_path, '--g', 'linewidth', 1); % Path of point B
plot(Cx_path, Cy_path, '--m', 'linewidth', 1); % Path of point C

% Hold the final configuration of the mechanism
plot([Ax Bx_path(end)], [Ay By_path(end)], 'r', 'LineWidth', 1);  % Final position of Link 2
plot([Bx_path(end) Cx_path(end)], [By_path(end) Cy_path(end)], 'b', 'LineWidth', 1);  % Final position of Link 3
plot(Bx_path(end), By_path(end), 'bo', 'MarkerFaceColor', 'w');   % Final Point B
plot(Cx_path(end), Cy_path(end), 'go', 'MarkerFaceColor', 'w');   % Final Point C

% Show the figure
hold off;