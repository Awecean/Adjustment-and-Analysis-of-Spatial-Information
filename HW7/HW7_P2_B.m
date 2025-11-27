% AASI HW07_P2 - Case II
% constraint h2,h3
% B11501037 Tao
%% Part 0 Clear
clear; clc; close all;
%% Part 1 Set parameter
%---------------------------
% Observed lengths
%---------------------------
s1 = 12.41;
s2 = 8.06;
s3 = 3.87;
s4 = 8.83;
s = [s1; s2; s3; s4];
%---------------------------
% Initial guess [x_D, x_E, y_D, y_E, a]
%---------------------------
X = [8; 12; 4; 8]; 
X_p = [1]; % slope a
X_all = [X;X_p];
%---------------------------
% Iteration parameters
%---------------------------
W = eye(4);
threshold = 1e-10;

disp("--- Start Iteration for Case II (Constraints h2, h3 and parameter a) ---");
%% Part 1 iterations
for iter = 1:20
    xD = X_all(1);  xE = X_all(2);
    yD = X_all(3);  yE = X_all(4);
    a  = X_all(5);
    
    % F(x)
    F = [
        (xE-0)^2 + (yE-4)^2 - s1^2;
        (xD-0)^2 + (yD-4)^2 - s2^2;
        (xD-8)^2 + (yD-0)^2 - s3^2;
        (xE-8)^2 + (yE-0)^2 - s4^2
    ];
    
    % B
    Bx = [
        0, 2*xE, 0, 2*yE-8;
        2*xD, 0, 2*yD-8, 0;
        2*xD-16, 0, 2*yD, 0;
        0, 2*xE-16, 0, 2*yE
    ];
    
    % Bn = [Bx | 0] 
    B = [Bx, zeros(4, 1)];

    % H = [h2,h3]^T
    H = [
        a*xE - yE - 4*a;
        a*xD - yD - 4*a
    ];
    
    % C = [D1 | D2]
    
    % D1 
    D1 = [
        0, a, 0, -1;  % (for h2)
        a, 0, -1, 0   % (for h3)
    ];
    
    % D2 
    D2 = [
        xE - 4; % dh2/da
        xD - 4  % dh3/da
    ];
    
    C = [D1, D2];
    
    % Construct the Augmented Normal Equation System
    
    N_mat = B' * W * B; 
    U_vec = -B' * W * F; 
    
    %left side 
    A_top = [N_mat, C'];
    A_bottom = [C, zeros(2, 2)];
    A = [A_top; A_bottom];
    %right side
    rhs = [U_vec; -H];
    
    % solve [dx, k]
    result = A \ rhs;
    
    dX = result(1:5); 
    % k = result(6:7); % Lagrangrian Multipliers
    
    %update
    X_all = X_all + dX;   % update
    
    % stop condition
    if norm(dX) < threshold
        disp("Converged!");
        break;
    end
end
%% Part 2 ouput
disp("--- Final Result for Case II (with parameter 'a') ---");
fprintf('Total Iterations: %d\n', iter);
disp("Solution [xD, xE, yD, yE, a] = ");
disp(X_all);
%% Part 3 ouput Cofactor matrix
N_inv = inv(N_mat(1:4,1:4)); % N^-1
Q_k = inv(D1 * N_inv * D1'); % %D1 is C in class [professor]
Q_XX = N_inv - N_inv * D1' * Q_k * D1 * N_inv; % Q_XX = N^-1 -N^-1 C^T M^_-1 C N^-1

disp("Cofactor Matrix Q_XX (for [xD, xE, yD, yE]) = 1e-3 *");
disp(Q_XX * 1e3);

%% Part 4 Check BDE on a line (h2, h3)
disp("--- Constraint Verification ---");
xD_final = X_all(1); xE_final = X_all(2); yD_final = X_all(3); yE_final = X_all(4); a_final = X_all(5);

h2_final_value = a_final * xE_final - yE_final - 4 * a_final;
h3_final_value = a_final * xD_final - yD_final - 4 * a_final;

fprintf('Final constraint value h2(X) = %.2e\n', h2_final_value);
fprintf('Final constraint value h3(X) = %.2e\n', h3_final_value);

tolerance = 1e-12; 
if all(abs([h2_final_value, h3_final_value]) < tolerance)
    disp('Conclusion: Both constraints h2 and h3 are satisfied. (Points D and E lie on the line y = a(x-4))');
else
    disp('Conclusion: The constraints are not satisfied within the tolerance.');
end
