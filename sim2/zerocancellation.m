clc; clear; 
%% 1. 真實值
a1_true = -1.6065;
a2_true =  0.6065;
b0_true =  0.1065;
b1_true =  0.0902;

%% 2.規格
m1 = -1.3205;
m2 =  0.4966;
Tm = 1 + m1 + m2;  

r1_true = b1_true / b0_true;
s0_true = (m1 - a1_true) / b0_true;
s1_true = (m2 - a2_true) / b0_true;
t0_true = Tm / b0_true;

%% 3. 模擬設定
N = 100;
time = 0:N-1;

uc = ones(1,N);
for i = 1:N/25
    if mod(i,2) == 0
        uc(1+(25*(i-1)):25*i) = -1;
    end
end

%% 4. RLS 初始化
lambda = 0.99;
P = diag([100, 100, 1, 1]);          
theta_hat = [0; 0; 0.01; 0.2];      
I = eye(4);

y1 = 0; y2 = 0;
u1 = 4; u2 = 0;

TH = zeros(4, N);  % theta
Y = zeros(1, N);   % y
U = zeros(1, N);   % u
RST = zeros(4,N);

TH(:,1) = theta_hat;
Y(1) = 0;      
U(1) = u1;      


%% 5. zero cancellation
for k = 2:N
    % 真實輸出
    y = -a1_true*y1 - a2_true*y2 + b0_true*u1 + b1_true*u2;

    % RLS 更新
    phi = [-y1; -y2; u1; u2];
    K = (P*phi) / (lambda + phi'*P*phi);
    e = y - phi'*theta_hat;
    theta_hat = theta_hat + K*e;
    P = (I - K*phi') * P / lambda;

    TH(:,k) = theta_hat;

    % 當前估值
    a1h = theta_hat(1);
    a2h = theta_hat(2);
    b0h = theta_hat(3);
    b1h = theta_hat(4);

    % RST
    r1 = b1h / b0h;
    s0 = (m1 - a1h) / b0h;
    s1 = (m2 - a2h) / b0h;
    t0 = Tm / b0h;

    % 記錄RST
    RST(:,k) = [r1; s0; s1; t0];

    % 控制律
    u = -r1*u1 + t0*uc(k) - s0*y - s1*y1;
    Y(k) = y;
    U(k) = u;

    y2 = y1;
    y1 = y;
    u2 = u1;
    u1 = u;
end

%% 6. 畫圖
figure('Name','EX3.4 Figure 3.4');
subplot(2,1,1)
plot(time, uc, 'k--', 'LineWidth', 1.2, 'DisplayName', '$u_c$'); hold on
plot(time, Y, 'b', 'LineWidth', 1.2, 'DisplayName', '$y$');
xlabel('Time');
ylabel(' $y$ and $u_c$', 'Interpreter', 'latex');
title('$u_c$ and $y$', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'best');
grid on;

subplot(2,1,2)
plot(time, U, 'k', 'LineWidth', 1.2);
xlabel('Time');
ylabel('Control $u$', 'Interpreter', 'latex');
title('Control Signal', 'Interpreter', 'latex');
grid on;

figure('Name','EX3.4 Figure 3.5');
subplot(2,1,1)
stairs(time, TH(1,:), 'b', 'LineWidth',1.2); hold on
stairs(time, TH(2,:), 'r', 'LineWidth',1.2);
yline(a1_true, 'b--', 'LineWidth', 1);
yline(a2_true, 'r--', 'LineWidth', 1);
xlabel('Time');
ylabel('Estimates');
title('Estimates of $a_1$ and $a_2$', 'Interpreter', 'latex');
legend({'$\hat{a}_1$', '$\hat{a}_2$', '$a_1\ true$', '$a_2\ true$'}, 'Interpreter', 'latex');
xlim([0 30]);
grid on;

subplot(2,1,2)
stairs(time, TH(3,:), 'b', 'LineWidth',1.2); hold on
stairs(time, TH(4,:), 'r', 'LineWidth',1.2);
yline(b0_true, 'b--', 'LineWidth', 1);
yline(b1_true, 'r--', 'LineWidth', 1);
xlabel('Time');
ylabel('Estimates');
title('Estimates of $b_0$ and $b_1$', 'Interpreter', 'latex');
legend({'$\hat{b}_0$', '$\hat{b}_1$', '$b_0\ true$', '$b_1\ true$'}, 'Interpreter', 'latex');
xlim([0 30]);
grid on;

figure('Name','RST');
hold on
yline(r1_true, '--r', 'LineWidth', 1.2);  
yline(s0_true, '--b', 'LineWidth', 1.2);  
yline(s1_true, '--g', 'LineWidth', 1.2); 
yline(t0_true, '--k', 'LineWidth', 1.2);
R1 = stairs(time, RST(1,:), 'r', 'LineWidth',1.2);
S0 = stairs(time, RST(2,:), 'b', 'LineWidth',1.2);
S1 = stairs(time, RST(3,:), 'g', 'LineWidth',1.2);
T0 = stairs(time, RST(4,:), 'k', 'LineWidth',1.2);
xlabel('Time');
ylabel('RST Estimates');
title('RST Parameter Estimates');
axis([0 20 -inf inf]);
legend([R1 S0 S1 T0], {'$r_1$', '$s_0$', '$s_1$', '$t_0$'}, 'Interpreter','latex','Location','best');
grid on;

%% 7. 平均值
avg_window = 20:N; 

a1_hat_avg = mean(TH(1, avg_window));
a2_hat_avg = mean(TH(2, avg_window));
b0_hat_avg = mean(TH(3, avg_window));
b1_hat_avg = mean(TH(4, avg_window));

r1_hat_avg = mean(RST(1, avg_window));
s0_hat_avg = mean(RST(2, avg_window));
s1_hat_avg = mean(RST(3, avg_window));
t0_hat_avg = mean(RST(4, avg_window));

fprintf('\nParameter Estimate Averages :\n');
fprintf('  a1_hat = %.6f   (true = %.6f)\n', a1_hat_avg, a1_true);
fprintf('  a2_hat = %.6f   (true = %.6f)\n', a2_hat_avg, a2_true);
fprintf('  b0_hat = %.6f   (true = %.6f)\n', b0_hat_avg, b0_true);
fprintf('  b1_hat = %.6f   (true = %.6f)\n', b1_hat_avg, b1_true);

fprintf('\nRST Parameter Averages :\n');
fprintf('  r1_hat = %.6f   (true = %.6f)\n', r1_hat_avg, r1_true);
fprintf('  s0_hat = %.6f   (true = %.6f)\n', s0_hat_avg, s0_true);
fprintf('  s1_hat = %.6f   (true = %.6f)\n', s1_hat_avg, s1_true);
fprintf('  t0_hat = %.6f   (true = %.6f)\n', t0_hat_avg, t0_true);