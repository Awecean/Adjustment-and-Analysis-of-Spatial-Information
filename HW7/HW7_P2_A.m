% AASI HW07_P2 - Case II
% constraint h1
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
% Initial guess [x_D, x_E, y_D, y_E]
%---------------------------
x = [8; 12; 4; 8];
%---------------------------
% Iteration parameters
%---------------------------
W = eye(4);
threshold = 1e-10;

disp("--- Start Iteration for Case II (Constrained) ---");
%% Part 1 iteration
for iter = 1:20
    xD = x(1);  xE = x(2);
    yD = x(3);  yE = x(4);
    
    %F(x)
    F = [
        (xE-0)^2 + (yE-4)^2 - s1^2;
        (xD-0)^2 + (yD-4)^2 - s2^2;
        (xD-8)^2 + (yD-0)^2 - s3^2;
        (xE-8)^2 + (yE-0)^2 - s4^2
    ];
    
    % B
    B = [
        0, 2*xE, 0, 2*yE-8;
        2*xD, 0, 2*yD-8, 0;
        2*xD-16, 0, 2*yD, 0;
        0, 2*xE-16, 0, 2*yE
    ];
    
    % H = h1
    h1 = -yD*xE - 4*yE + 4*yD + yE*xD;
    h = h1; 

    % C = D1
    D1 = [yE, -yD, 4-xE, xD-4];
    
    % Construct the Augmented Normal Equation System 
    
    % M = N_mat
    
    N_mat = B' * W * B;
    U_vec = -B' * W * F;
    
    %left side
    A_top = [N_mat, D1'];
    A_bottom = [D1, 0];
    A = [A_top; A_bottom];
    %right side
    rhs = [U_vec; -h];
    
    % solve [dx, k]
    result = A \ rhs;
    
    dx = result(1:4);
    k = result(5); % Lagrangrian Multiplier
    
    %update
    x = x + dx;   % update
    
    % stop condition
    if norm(dx) < threshold
        disp("Converged!");
        break;
    end
end
%% Part 2 ouput x
disp("--- Final Result for Case II ---");
disp("Solution [xD, xE, yD, yE] = ");
disp(x);

%% Part 3 ouput Cofactor matrix
N_inv = inv(N_mat);
Q_k = inv(D1 * N_inv * D1'); % Q_k is 1x1
Q_xx = N_inv - N_inv * D1' * Q_k * D1 * N_inv;

disp("Cofactor Matrix Q_xx = 1e-3 * ");
disp(Q_xx * 1e3);


%% Part 4 Check BDE on a line (h1)
disp("--- Constraint Verification ---");

xD_final = x(1);
xE_final = x(2);
yD_final = x(3);
yE_final = x(4);

h1_final_value = -yD_final * xE_final - 4 * yE_final + 4 * yD_final + yE_final * xD_final;

fprintf('Final constraint value h1(x) = %.2e\n', h1_final_value);

tolerance = 1e-12; 

if abs(h1_final_value) < tolerance
    disp('Conclusion: The solution [D, E] is confirmed to be on a straight line with point B(4, 0). (Constraint satisfied)');
else
    disp('Conclusion: The constraint h1 is not satisfied within the tolerance.');
end
% due to there are multiply of parameter, the result h1 might not equal to
% 0, I thought.