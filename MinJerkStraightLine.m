% Modelling of reaching movements using minimum jerk trajectories of human
% reaching movements.
% The different equations are from Flash & Hogan, 1986
% @Antoine de Comite - 05/04/2021

%% 0. Hyperparameters for the straight line movements

clear all; clc; close all;
% Initial position
x0 = 0; y0 = 0;
% Final position
xf=0.25; yf=0.25;
% Time vector 
t_vector = linespace(0,1,1000);
tau = t_vector/t_vector(end);

%% 1. Modelling of straight trajectories

% Computation of the position 
x_position = zeros(length(t_vector),1);
y_position = zeros(length(t_vector),1);

for ii = 1 : length(t_vector)
    x_position(ii) = x0 + (xf-x0)*(6*tau(ii)^5-15*tau(ii)^4+10*tau(ii)^3);
    y_position(ii) = y0 + (yf-y0)*(6*tau(ii)^5-15*tau(ii)^4+10*tau(ii)^3);
end

% Computation of hte speed 
x_speed = zeros(length(x_position),1); y_speed = zeros(length(y_position),1);
x_speed(3:end-2) = (-x_position(5:end)+8*x_position(4:end-1)-8*x_position(2:end-3)+x_position(1:end-4))/(12*(t_vector(2)-t_vector(1)));
y_speed(3:end-2) = (-y_position(5:end)+8*y_position(4:end-1)-8*y_position(2:end-3)+y_position(1:end-4))/(12*(t_vector(2)-t_vector(1)));

figure('Name','Minimum jerk, straight line unperturbed'); hold on; 
subplot(1,2,1); hold on; title('Position');
set(gca,'Color','none');
plot(x_position,y_position,'b','LineWidth',2);
plot(x0,y0,'g.','MarkerSize',25);
plot(xf,yf,'r.','MarkerSize',25);
xlabel('x-position'); ylabel('y-position');

subplot(1,2,2); hold on; title('x- and y-speeds');
set(gca,'Color','none');
plot(tau,x_speed,'b','LineWidth',2);
plot(tau,y_speed,'b:','LineWidth',2);
legend('x-speed','y-speed'); xlabel('Relative time [t/tf]'); ylabel('Movement speed [m/s]');

