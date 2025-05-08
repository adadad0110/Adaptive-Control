clc; clear;

%% 1.  真實值
a1_true = -1.6065;
a2_true =  0.6065;
b0_true =  0.1065;
b1_true =  0.0902;

%% 2.  規格
m1 = -1.3205;
m2 =  0.4966;
Tm = 1 + m1 + m2;

den_true = b1_true^2 - a1_true*b0_true*b1_true + a2_true*b0_true^2;
r1_true  = ( b1_true^2*(m1 - a1_true)- b0_true*b1_true*(m2 - a2_true) ) / den_true;
s0_true  = ( m1 - a1_true - r1_true ) / b0_true;
s1_true  = -a2_true * r1_true / b1_true;
t0_true  = (1 + m1 + m2) / (b0_true + b1_true);
%% 3. 模擬設定
N    = 800;
time = 0:N-1;

uc = ones(1, N);
for i = 1:(N/25)
    if mod(i,2) == 0
        uc((25*(i-1)+1):(25*i)) = -1;
    end
end

%% 4. RLS 初始化
lambda      = 0.99;
P           = diag([100,100,1,1]);
theta_hat   = [0;0;0.01;0.2];
I4          = eye(4);

y1 = 0; y2 = 0;
u1 = 1; u2 = 0;

TH  = zeros(4, N);
Y   = zeros(1, N);
U   = zeros(1, N);
RST = zeros(4, N);

TH(:,1) = theta_hat;
Y(1)    = 0;
U(1)    = u1;

%% 5. no zero cancellation
for k = 2:N
    % 真實輸出
    y = -a1_true*y1 - a2_true*y2 + b0_true*u1 + b1_true*u2;

    % RLS 更新
    phi = [-y1; -y2; u1; u2];
    K = (P*phi) / (lambda + phi'*P*phi);
    e = y - phi'*theta_hat;
    theta_hat = theta_hat + K*e;
    P = (I4 - K*phi')*P / lambda;

    TH(:,k) = theta_hat;

    % 當前估值
    a1h = theta_hat(1);
    a2h = theta_hat(2);
    b0h = theta_hat(3);
    b1h = theta_hat(4);

    % no zero cancellation RST
    den = b1h^2 - a1h*b0h*b1h + a2h*b0h^2;
    r1  = ( b1h^2*(m1 - a1h) - b0h*b1h*(m2 - a2h) ) / den;
    s0  = ( m1 - a1h - r1 ) / b0h;
    s1  = -a2h * r1 / b1h;
    t0  = (1 + m1 + m2) / (b0h + b1h);

    RST(:,k) = [r1; s0; s1; t0];

    % 控制律
    u       = -r1*u1 + t0*uc(k) - s0*y - s1*y1;
    Y(k)    = y;
    U(k)    = u;

    y2 = y1;  
    y1 = y;
    u2 = u1;  
    u1 = u;
end

%% 6. 畫圖
figure('Name','EX3.5 Figure 3.6');
subplot(2,1,1);
plot(time, uc,'k--','LineWidth',1.2); hold on;
plot(time, Y, 'b-','LineWidth',1.2); grid on;
legend('$u_c$','$y$','Interpreter','latex');
title('Process Output vs Command','Interpreter','latex');
xlim([0 100]);

subplot(2,1,2);
plot(time, U,'k-','LineWidth',1.2); grid on;
title('Control Signal','Interpreter','latex');
xlim([0 100]);

figure('Name','EX3.5 Figure 3.7');
subplot(2,1,1);
stairs(time, TH(1,:), 'b-'); hold on;
stairs(time, TH(2,:), 'r-');
yline(a1_true,'b--'); yline(a2_true,'r--');
legend('$\hat a_1$','$\hat a_2$','$a_1$ true','$a_2$ true','Interpreter','latex');
title('Estimates of $a_i$','Interpreter','latex'); grid on;

subplot(2,1,2);
stairs(time, TH(3,:), 'b-'); hold on;
stairs(time, TH(4,:), 'r-');
yline(b0_true,'b--'); yline(b1_true,'r--');
legend('$\hat b_0$','$\hat b_1$','$b_0$ true','$b_1$ true','Interpreter','latex');
title('Estimates of $b_j$','Interpreter','latex'); grid on;

figure('Name','RST'); 
hold on;
yline(r1_true, '--r', 'LineWidth', 1.2);  
yline(s0_true, '--b', 'LineWidth', 1.2);  
yline(s1_true, '--g', 'LineWidth', 1.2); 
yline(t0_true, '--k', 'LineWidth', 1.2);
stairs(time, RST(1,:),'r-'); 
stairs(time, RST(2,:),'b-');
stairs(time, RST(3,:),'g-'); 
stairs(time, RST(4,:),'k-');
legend('$r_1$','$s_0$','$s_1$','$t_0$','Interpreter','latex');
title('RST Parameters (no zero cancellation)','Interpreter','latex');
grid on;

%% 7. 平均值
avg_window = 400:800; 

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