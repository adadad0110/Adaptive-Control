clc; clear; close all;

%% 可調整參數
order = 4;             % 4, 6, 8
alpha = 10;             % a>0
r = 0.5;               % 2>r>0
x_raw = linspace(0, 3, 2000)';   
x = x_raw / max(x_raw);         
N = length(x);
u = x;
y = exp(-x_raw);               

%% 特徵矩陣 Phi
Phi = ones(N, order + 1);
for i = 1:order
    Phi(:, i+1) = u.^i;
end

%% 初始 Projection
theta_hat = zeros(order + 1, 1);
theta_hat_history = zeros(N, order + 1);
y_hat = zeros(N, 1);

%% Projection Approximation
for t = 1:N
    phi_t = Phi(t, :)';
    y_t = y(t);
    e_t = y_t - phi_t' * theta_hat;
    ccc = alpha + phi_t' * phi_t;
    theta_hat = theta_hat + r * phi_t * e_t / ccc;
    y_hat(t) = phi_t' * theta_hat;
    theta_hat_history(t, :) = theta_hat';
end

theta_true = zeros(order+1, 1);
for i = 0:order
    theta_true(i+1) = ((-1)^i) / factorial(i);
end

%% 圖1 真實 估計輸出
figure;
subplot(1,2,1);
plot(x_raw, y, 'b', 'LineWidth', 2, 'DisplayName', '$y = \exp(-x)$'); hold on;
plot(x_raw, y_hat, 'r--', 'LineWidth', 2, ...
    'DisplayName', ['$\hat{y}$ (Projection, order = ', num2str(order), ')']);
xlabel('x'); ylabel('y');
title(['Output Comparison (Projection, Order = ', num2str(order), ')'], 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on;

%% 圖2 theta
subplot(1,2,2);
plot(1:N, theta_hat_history, 'LineWidth', 1.5); hold on;
legend(arrayfun(@(i) ['$\theta_', num2str(i-1), '$'], 1:(order+1), 'UniformOutput', false), ...
       'Interpreter', 'latex');
xlabel('Iteration'); ylabel('$\hat{\theta}$', 'Interpreter', 'latex');
title('Estimated Parameters $\hat{\theta}(t)$ over Time', 'Interpreter', 'latex');
grid on;

%% 圖3 true & estimated
theta_true = zeros(order+1, 1);
for i = 0:order
    theta_true(i+1) = ((-1)^i) / factorial(i);
end

figure;
rows = ceil(sqrt(order + 1));
cols = ceil((order + 1) / rows);

for i = 1:(order+1)
    subplot(rows, cols, i);

    % 估算
    plot(1:N, theta_hat_history(:, i), 'b', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Estimated $\\hat{\\theta}_{%d}$', i-1));
    hold on;

    % 真實
    yline(theta_true(i), 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('True $\\theta_{%d}$ = %.4f', i-1, theta_true(i)), ...
        'Interpreter', 'latex');

    xlabel('Iteration');
    ylabel(['$\theta_', num2str(i-1), '(t)$'], 'Interpreter', 'latex');
    title(['Parameter $\theta_', num2str(i-1), '$ over Time (Projection)'], 'Interpreter', 'latex');

    legend('Interpreter', 'latex', 'Location', 'best');
    grid on;
end
