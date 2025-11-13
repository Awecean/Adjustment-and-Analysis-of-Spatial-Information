% AASI HW06
% B11501036 Tao

%% Part 0
clear; clc; close all;

%% Part 1 parameter setting
x_L = [8.928, 14.808, 14.719, 18.526, 19.230];
y_L = [3.739, 16.608, 6.759, 17.255, 10.057];
x_U = [1.477, 2.107, 4.010, 4.245, 5.607];
y_U = [1.229, 4.109, 1.562, 3.636, 1.615];

%% Part 2 - Problem 1.a
B = [x_U',y_U',ones(5,1),zeros(5,3);
    zeros(5,3),x_U',y_U',ones(5,1)];
f = [x_L,y_L]';
l = f;
Q_ff = 10^(-4)*eye(10);
w = inv(Q_ff);
%w = *diag(ones(1,10));
K_hat = inv(B' * w * B)*(B' * w * f);
v = f-B*K_hat;
l_hat = l-v;
disp('--------------Problem 1--------------')
disp('K_hat^T');disp(K_hat');
disp('v^T');disp(v');
disp('l_hat^T');disp(l_hat');
Q_KK = inv(B' * w * B);
disp('Q_{hat(K)hat(K)}'); disp(Q_KK);
disp('-------------------------------------')

%% Part 3 - Problem 1.b

K0 = K_hat;
l0 = l_hat;

W = [1e4*eye(5),zeros(5,5); %x-y
    zeros(5,5), 1e2*(eye(5))]; %X-Y

deltaK = 5*ones(6,1);    % set initial difference step
iter = 0;
tol = 1e-8;

while norm(deltaK) > tol 
    iter = iter + 1;

    X = x_U(:);
    Y = y_U(:);
    x = l0(1:5);
    y = l0(6:10);

    % A and B
    A = [-eye(5), zeros(5), K0(1)*eye(5), K0(2)*eye(5);
         zeros(5), -eye(5), K0(4)*eye(5), K0(5)*eye(5)];

    B = [X, Y, ones(5,1), zeros(5,3);
         zeros(5,3), X, Y, ones(5,1)];

    % f
    f_val = [K0(1)*X + K0(2)*Y + K0(3) - x;
             K0(4)*X + K0(5)*Y + K0(6) - y];

    % \Delta K, v
    deltaK = -(B' * W * B) \ (B' * W * f_val);
    v = -A * ((A' * W * A) \ (A' * W * (f_val - B * deltaK)));

    % Update
    K0 = K0 + deltaK;
    l0 = l0 + v;

    fprintf('Iteration %d, ||deltaK|| = %.3e\n', iter, norm(deltaK));
end

K_hat2 = K0;
l_hat2 = l0;

% Covariance matrix

Q_KK_2 = inv(B'*W*B);

disp('--------------Problem 2--------------')
disp('Final estimated parameters K_hat^T = ')
disp(K_hat2')
disp('Covariance matrix 1e3*Q_{KK} = ') % to compare easily
disp(Q_KK_2*1e3)
disp('-------------------------------------')

%% Part 4 Problem 2 comparison

sqrt_diag_1 = sqrt(diag(Q_KK));
disp(sqrt_diag_1');
sqrt_diag_2 = sqrt(diag(Q_KK_2));
disp(sqrt_diag_2');
% the parameter e and f 's uncertainty is greater in second case.
