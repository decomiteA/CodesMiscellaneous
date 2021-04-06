% Modelling of via point movements with minimum jerk.
% Equations from Flash & Hogan 1986
%
% @Antoine de Comite

close all;
% Some symbolic variable definition
syms tau_1 tf x0 x1 xf y0 y1 yf float
tau_1 = sym('tau_1','real');
tf = sym('tf','real');
x0 = sym('x0','real');
x1 = sym('x1','real');
xf = sym('xf','real');
y0 = sym('y0','real');
y1 = sym('y1','real');
yf = sym('yf','real');

% Set here the values that you want for the different parameters
% Initial position
x0 = 0; y0 = 0;
% Final position
xf = 0; yf = 0.3;
% Via point position
x1 = 0.1; y1 = 0.25;
% Time vector
t_vector = linspace(0,1,1000);
tau = t_vector/t_vector(end);
tf = t_vector(end);
% Solving the polynomial as described in Hogan & Flash and Sharkawy's
% papers.


A1x = (xf-x0)*(300*tau_1^5-1200*tau_1^4+1600*tau_1^3) + tau_1^2*(-720*xf+120*x1+600*x0) + (x0-x1)*(300*tau_1-200);
A1y = (yf-y0)*(300*tau_1^5-1200*tau_1^4+1600*tau_1^3) + tau_1^2*(-720*yf+120*y1+600*y0) + (y0-y1)*(300*tau_1-200);
A2x = (xf-x0)*(120*tau_1^5-300*tau_1^4+200*tau_1^3) - 20*(x1-x0);
A2y = (yf-y0)*(120*tau_1^5-300*tau_1^4+200*tau_1^3) - 20*(y1-y0);
A3 = 60*tau_1^7 - 210*tau_1^6 + 240*tau_1^5 - 90*tau_1^4;
A4 = 60*tau_1^3 - 30*tau_1^2 - 30*tau_1^4;
FonctionToSolve = A3*A2x*A2x + A3*A2y*A2y + tau_1^3 * (A1x*A2x*A4) + tau_1^3 * (A1y*A2y*A4);
solvet = vpasolve(FonctionToSolve,tau_1);


% The following might not be the more robust... Tbh I do not know how
% robust this is and it might need some tweak... I had to select the only
% real-valued root of the polynomial FonctionToSolve.

tau_1 = real(solvet(6)); tau_1 = double(tau_1);
[~,idx_tau] = min(abs(tau-tau_1));

% Let us now define some constants and cast them in the right format

% Constant definition
c1 = 1/(tf^5*tau_1^2*(1-tau_1)^5) * ((xf-x0)*(300*tau_1^5-1200*tau_1^4+1600*tau_1^3)+tau_1^2*(-720*xf+120*x1+600*x0)+(x0-x1)*(300*tau_1-200));
c2 = 1/(tf^5*tau_1^2*(1-tau_1)^5) * ((yf-y0)*(300*tau_1^5-1200*tau_1^4+1600*tau_1^3)+tau_1^2*(-720*yf+120*y1+600*y0)+(y0-y1)*(300*tau_1-200));
pi_1 = 1/(tf^5*tau_1^5*(1-tau_1)^5) * ((xf-x0)*(120*tau_1^5-300*tau_1^4+200*tau_1^3)-20*(x1-x0));
pi_2 = 1/(tf^5*tau_1^5*(1-tau_1)^5) * ((yf-y0)*(120*tau_1^5-300*tau_1^4+200*tau_1^3)-20*(y1-y0));

% Constant casting (not 100% necessary but allows to circumvent some
% errors)


pi_1 = double(pi_1); pi_2 = double(pi_2);
c1 = double(c1); c2 = double(c2);


% Computation of the x- and y-positions with respect to time for both
% segments

% 1st segment before via point
x_before = zeros(length(1:idx_tau),1);
y_before = zeros(length(1:idx_tau),1);
for ii = 1 : length(x_before)
   x_before(ii) = (tf)^5/720*(pi_1*(tau_1^4*(15*tau(ii)^4-30*tau(ii)^3)+tau_1^3*(80*tau(ii)^3-30*tau(ii)^4)-60*tau(ii)^3*tau_1^2+30*tau(ii)^4*tau_1-6*tau(ii)^5)+c1*(15*tau(ii)^4-10*tau(ii)^3-6*tau(ii)^5))+x0; 
   y_before(ii) = (tf)^5/720*(pi_2*(tau_1^4*(15*tau(ii)^4-30*tau(ii)^3)+tau_1^3*(80*tau(ii)^3-30*tau(ii)^4)-60*tau(ii)^3*tau_1^2+30*tau(ii)^4*tau_1-6*tau(ii)^5)+c2*(15*tau(ii)^4-10*tau(ii)^3-6*tau(ii)^5))+y0; 
end

% 2nd segment after via point
x_after = zeros(length(idx_tau+1:length(tau)),1);
y_after = zeros(length(idx_tau+1:length(tau)),1);
for ii = 1 : length(x_after)
   iii = ii+idx_tau;
   x_after(ii) = (tf)^5/720*(pi_1*(tau_1^4*(15*tau(iii)^4-30*tau(iii)^3)+tau_1^3*(80*tau(iii)^3-30*tau(iii)^4)-60*tau(iii)^3*tau_1^2+30*tau(iii)^4*tau_1-6*tau(iii)^5)+c1*(15*tau(iii)^4-10*tau(iii)^3-6*tau(iii)^5))+x0+pi_1*((tf^5*(tau(iii)-tau_1)^5)/(120));
   y_after(ii) = (tf)^5/720*(pi_2*(tau_1^4*(15*tau(iii)^4-30*tau(iii)^3)+tau_1^3*(80*tau(iii)^3-30*tau(iii)^4)-60*tau(iii)^3*tau_1^2+30*tau(iii)^4*tau_1-6*tau(iii)^5)+c2*(15*tau(iii)^4-10*tau(iii)^3-6*tau(iii)^5))+y0+pi_2*((tf^5*(tau(iii)-tau_1)^5)/(120));
end


x_total = [x_before;x_after];
y_total = [y_before;y_after];
x_speed = zeros(length(x_total),1); y_speed = zeros(length(y_total),1);
x_speed(3:end-2) = (-x_total(5:end)+8*x_total(4:end-1)-8*x_total(2:end-3)+x_total(1:end-4))/(12*(t_vector(2)-t_vector(1)));
y_speed(3:end-2) = (-y_total(5:end)+8*y_total(4:end-1)-8*y_total(2:end-3)+y_total(1:end-4))/(12*(t_vector(2)-t_vector(1)));

% Graphical representation of the movement 
figure('Name','Representation for the curved movements'); hold on;
subplot(1,2,1); hold on; set(gca,'Color','none'); title('Position');
plot(x_total,y_total,'r','LineWidth',3);
plot(x0,y0,'g.','MarkerSize',20);
plot(x1,y1,'k.','MarkerSize',20);
plot(xf,yf,'r.','MarkerSize',20);

subplot(1,2,2); hold on; set(gca,'Color','none'); title('Speed');
plot(tau,x_speed,'b','Linewidth',3);
plot(tau,y_speed,'b:','LineWidth',3);
xlabel('Relative time [tau/tf]'); ylabel('Speed [m/s]');