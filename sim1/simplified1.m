clear; clc; close;

%% 可調整參數
order = 6;                  % 4, 6, 8
N=2000;
x = linspace(0, 3, N)';  
y = exp(-x);                
delta = 0.01;                
alpha = 10;                 

%% 初始
phi = zeros(N, order);      
for i = 1:order
    phi(:, i) = x.^i;
end
phi = [ones(N,1), phi];     
n_params = size(phi, 2);

theta_hat = zeros(n_params, N);   
y_hat = zeros(N,1);               
P = alpha * eye(n_params);        
theta = zeros(n_params, 1);       

%% Simplified RLS
for t = 1:N
    phi_t = phi(t, :)';                           
    y_t = y(t);                                   
    K = P * phi_t / (delta + phi_t' * P * phi_t); 
    e_t = y_t - phi_t' * theta;                   
    theta = theta + K * e_t;                      
    P = (eye(n_params) - K * phi_t') * P;         
    theta_hat(:, t) = theta;                     
    y_hat(t) = phi_t' * theta;                    
end

%% 圖1 真實 估計輸出
figure;
subplot(1,2,1);
plot(x, y, 'b', 'LineWidth', 2, 'DisplayName', '$y = \exp(-x)$'); hold on;
plot(x, y_hat, 'r--', 'LineWidth', 2, 'DisplayName', ['$\hat{y}$ (order = ', num2str(order), ')']);
xlabel('x');
ylabel('y');
title(['Output Comparison (Order = ', num2str(order), ')'], 'Interpreter', 'latex');
legend('Interpreter', 'latex');
grid on;

%% 圖2 theta
subplot(1,2,2);
plot(1:N, theta_hat', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('$\hat{\theta}$', 'Interpreter', 'latex');
title('Estimated Parameters $\hat{\theta}(t)$ over Time', 'Interpreter', 'latex');
legend(arrayfun(@(i) ['$\theta_', num2str(i-1), '$'], 1:(order+1), 'UniformOutput', false), ...
       'Interpreter', 'latex', 'Location', 'best');
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
    plot(1:N, theta_hat(i, :), 'b', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Estimated $\\hat{\\theta}_{%d}$', i-1));
    hold on;

    % 真實
    yline(theta_true(i), 'r--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('True $\\theta_{%d}$ = %.4f', i-1, theta_true(i)), ...
        'Interpreter', 'latex');

    xlabel('Iteration');
    ylabel(['$\theta_', num2str(i-1), '(t)$'], 'Interpreter', 'latex');
    title(['Parameter $\theta_', num2str(i-1), '$ over Time'], 'Interpreter', 'latex');

    legend('Interpreter', 'latex', 'Location', 'best');
    grid on;
end
