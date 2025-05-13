clc; clear; close all;

%% parameters
a = 1;
b = 0.5;
am = 2;
bm = 2;
theta1_true = bm / b;      % = 4
theta2_true = (am - a) / b; % = 2
gamma = 5;  

%% settings
T  = 1000;    
dt = 0.1;    
time = 0:dt:T;
n = length(time);

%% Input
uc = ones(1,n);
uc(mod(time,20) >= 10) = -1;

y       = zeros(1,n);
ym      = zeros(1,n);  
u       = zeros(1,n);   
theta1  = zeros(1,n);  
theta2  = zeros(1,n);  

%% loop
for k = 2:n
    dym = -am * ym(k-1) + bm * uc(k-1);
    ym(k) = ym(k-1) + dt * dym;

    % Control law
    u(k) = theta1(k-1) * uc(k-1) - theta2(k-1) * y(k-1);

    % Plant update
    dy = -a * y(k-1) + b * u(k);
    y(k) = y(k-1) + dt * dy;

    e = y(k) - ym(k);

    F    = 1 / (1 + dt * am);
    phi1 = F * uc(k-1);
    phi2 = F * y(k-1);

    theta1(k) = theta1(k-1) - dt * gamma * phi1 * e;
    theta2(k) = theta2(k-1) + dt * gamma * phi2 * e;
end

%% Figure 5.5
figure(1);
subplot(2,1,1);
plot(time, ym, 'b', time, y, 'r', 'LineWidth', 1.2);
legend('y_m', 'y');
title(sprintf('Figure 5.5 (\\gamma=%.2f)', gamma));
axis([0, 100, -1.5, 1.5]);
grid on;

subplot(2,1,2);
plot(time, u, 'k', 'LineWidth', 1.2);
axis([0, 100, -5.5, 5.5]);
grid on;

%% Figure 5.6
figure(2);
subplot(2,1,1);
plot(time, theta1, 'LineWidth', 1.2); hold on;
plot([0 T], [theta1_true theta1_true], 'k--', 'LineWidth', 1);
legend('Estimate', 'True');
title(sprintf('Figure 5.6: \\theta_1 (\\gamma=%.2f)', gamma));
axis([0, 100, 0, 4.3]);
grid on;

subplot(2,1,2);
plot(time, theta2, 'LineWidth', 1.2); hold on;
plot([0 T], [theta2_true theta2_true], 'k--', 'LineWidth', 1);
legend('Estimate', 'True');
title(sprintf('Figure 5.6: \\theta_2 (\\gamma=%.2f)', gamma));
axis([0, 100, -1, 2.3]);
grid on;

%% Figure 5.7
figure(3);
plot(theta1, theta2, 'LineWidth', 1.2); hold on;
x_ref = [1, theta1_true];
y_ref = x_ref - (a/b);
plot(x_ref, y_ref, 'k--', 'LineWidth', 1);
%convergence point
plot(theta1_true, theta2_true, 'ko', 'MarkerFaceColor', 'k');
axis([0, 4.5, -1, 2.5]);
legend(' ', 'θ_2 = θ_1 - a/b', 'Convergence point');
title(sprintf('Relation θ_2 & θ_1 (\\gamma=%.1f)', gamma));
grid on;
