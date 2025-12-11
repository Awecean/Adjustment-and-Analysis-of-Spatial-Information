% AASI HW08 P1 - Case II: Constraint (h1)
%% Part 0: Initialization
clear; clc;
%% Part 1 parameter declare
xA = 0; yA = 4; xB = 4; yB = 0; xC = 8; yC = 0; 
s = [12.41; 8.06; 3.87; 8.83]; % s1, s2, s3, s4

% Uncertainty
sigma_s = 0.02;   %distance    
sigma_c = 0.01;   %coordinate

x = [8; 12; 4; 8]; 
threshold = 1e-10;
max_iter = 20;

Q_LL_s = (sigma_s^2) * eye(4);      
     
Q_XX = 1e8 * eye(4);

% dimension of matrix
n_L = 4; 
n_c = 1; 
n_x = 4; 

A = diag(-2*s);             
Ac = eye(n_c);    

%% Part 2 Iteration
disp("--- Start Iteration---");
for iter = 1:max_iter
    xD = x(1);  xE = x(2); yD = x(3);  yE = x(4);

    % B = dF/dX
    B = [
        0, 2*(xE-xA), 0, 2*(yE-yA);
        2*(xD-xA), 0, 2*(yD-yA), 0;
        2*(xD-xC), 0, 2*(yD-yC), 0;
        0, 2*(xE-xC), 0, 2*(yE-yC)
    ];
    
    % F: f1 to f4
    F = [
        (xE-xA)^2 + (yE-yA)^2 - s(1)^2;
        (xD-xA)^2 + (yD-yA)^2 - s(2)^2;
        (xD-xC)^2 + (yD-yC)^2 - s(3)^2;
        (xE-xC)^2 + (yE-yC)^2 - s(4)^2
    ];
    % J_C = dF/dcoord(C)
    J_C = [
    0, 0;
    0, 0;
    -2*(xD-xC), -2*(yD-yC);
    -2*(xE-xC), -2*(yE-yC)
    ];
    Q_LL_c= J_C *(sigma_c^2) * eye(2)* J_C';
    Q_LL = Q_LL_s + Q_LL_c;
    % H1 
    H1 = yE * (xD - xB) - yD * (xE - xB); 
    
    % C = dH1/dx 
    C = [yE, -(yD), (xB - xE), (xD - xB)]; 
    
    %
    A_T = [A, zeros(n_L, n_c), zeros(n_L, n_x);
           zeros(n_c, n_L), Ac, zeros(n_c, n_x);
           zeros(n_x, n_L), zeros(n_x, n_c), eye(n_x)]; 
    
    
    B_T = [B; C; eye(n_x)]; % B,C,I
    
    f_T = [F; H1; zeros(n_x, 1)]; 
    
    Q_c = sigma_c^2*((yD-yE)^2+(xE-xD)^2); % as mention in report.
    
    Q_T = blkdiag(Q_LL, Q_c, Q_XX); 
    
    Q_eT = A_T * Q_T * A_T';
    W_eT = pinv(Q_eT);

    % N_T = B_T^T * W_eT * B_T
    N_T = B_T' * W_eT * B_T;
    
    % t_T = B_T^T * W_eT * f_T
    t_T = B_T' * W_eT * f_T;
    
    % Solve for Delta (dx)
    dx = pinv(N_T) * t_T;

        
    % Update
    x_old = x;
    x = x - dx;   
    
    % Check
    if norm(dx) < threshold
        disp("Converged!");
        break;
    end
end
%% Part 3 print hat(x)
disp("--- Final Result for Case II---");
fprintf('Total Iterations: %d\n', iter);
disp("Estimated Coordinates x_hat [xD, xE, yD, yE] = ");
disp(x);

%% Part 4. Compute Q_xx
Sigma_xx = inv(N_T); 
disp("Final Covariance Matrix Sigma_xx (Q_xx) = 1e-3 *");
disp(Sigma_xx * 1e3);
%% for comparison
disp(norm(Sigma_xx));