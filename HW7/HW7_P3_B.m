% AASI HW07_P2 - Case V: Constrained Adjustment (h1, h6, h7 and parameter r)
%% Part 0 Clear
clear; clc;
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
% Initial guess [xD, xE, yD, yE, r]
r0 = sqrt((8-4)^2 + (4-0)^2); %guess r by BD(D is guessed)
X = [8; 12; 4; 8; r0]; % X = [xD; xE; yD; yE; r]
% Iteration parameters
W = eye(4);
threshold = 1e-10;

disp("--- Start Iteration for Case V (h1, h6, h7, parameter r) ---");
%% Part 1 iteration
for iter = 1:20
    xD = X(1);  xE = X(2);
    yD = X(3);  yE = X(4);
    r  = X(5);
    
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
    B = [Bx, zeros(4, 1)]; % Bn = [Bx | 0]

    % H = [h1;h6;h7]
    h1 = -yD*xE - 4*yE + 4*yD + yE*xD;
    h6 = (yE - yD)^2 + (xE - xD)^2 - r^2;
    h7 = (0 - yD)^2 + (4 - xD)^2 - r^2;
    H = [h1; h6; h7];
    
    %C = [D1, D2]
    

    D1= [yE, -yD, 4-xE, xD-4;
    2*(xD-xE), 2*(xE-xD), 2*(yD-yE), 2*(yE-yD);
    2*(xD-4), 0, 2*yD, 0];

    D2 = [0; -2*r; -2*r];
    
    C = [D1, D2]; 
    
    %Construct the Augmented Normal Equation System
    
    N_mat = B' * W * B; 
    U_vec = -B' * W * F; 
    
    %left side 
    A_top = [N_mat, C'];
    A_bottom = [C, zeros(3, 3)];
    A = [A_top; A_bottom]; % A is 8x8
    %right side
    rhs = [U_vec; -H]; % rhs is 8x1
    
    %solve [dx,k]
    result = A \ rhs;
    
    dX = result(1:5); 
    
    %update
    X = X + dX;
    
    % stop condition
    if norm(dX) < threshold
        disp("Converged!");
        break;
    end
end
%% Part 2 ouput x
disp("--- Final Result for Case V (h1, h6, h7, parameter r) ---");
fprintf('Total Iterations: %d\n', iter);
disp("Solution [xD, xE, yD, yE, r] = ");
disp(X);

%% Part 3 ouput Cofactor matrix
A_inv = inv(A); 
Q_XX = A_inv(1:5, 1:5);
Q_xx = A_inv(1:4, 1:4); 

disp("Cofactor Matrix Q_xx (for [xD, xE, yD, yE]) = 1e-3 *");
disp(Q_xx * 1e3);

%% Part 4 Check constraints (h1, h6, h7)
disp("--- Constraint Verification ---");
xD_final = X(1); xE_final = X(2); yD_final = X(3); yE_final = X(4); r_final = X(5);

h1_final = -yD_final*xE_final - 4*yE_final + 4*yD_final + yE_final*xD_final;
h6_final = (yE_final - yD_final)^2 + (xE_final - xD_final)^2 - r_final^2;
h7_final = (0 - yD_final)^2 + (4 - xD_final)^2 - r_final^2;

fprintf('Final constraint value h1(X) = %.2e\n', h1_final);
fprintf('Final constraint value h6(X) = %.2e\n', h6_final);
fprintf('Final constraint value h7(X) = %.2e\n', h7_final);

tolerance = 1e-12; 
if all(abs([h1_final, h6_final, h7_final]) < tolerance)
    disp('Conclusion: All constraints are satisfied.');
else
    disp('Conclusion: The constraints are NOT satisfied within the tolerance.');
end