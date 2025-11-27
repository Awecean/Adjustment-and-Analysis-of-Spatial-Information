% AASI HW07_P3
% constraints: h1, h4
% B11501037 Tao
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
% Initial guess [x_D, x_E, y_D, y_E]
%---------------------------
x = [8; 12; 4; 8]; % x = [xD; xE; yD; yE]
%---------------------------
% Iteration parameters
%---------------------------
W = eye(4);
threshold = 1e-10;

disp("--- Start Iteration for Case IV (Constraints h1, h4 ONLY) ---");
%% Part 1 iterations
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
    
    % H = [h1;h4]
    H = [
        -yD*xE - 4*yE + 4*yD + yE*xD; % h1
        2*xD - xE - 4                 % h4 
    ]; 
    
    %C = [D1; D2]
    C = [
        yE, -yD, 4-xE, xD-4; % dh1/dx
        2, -1, 0, 0          % dh4/dx
    ]; 
    
    %Construct the Augmented Normal Equation System 
    
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
    dx = result(1:4);
    
    %update
    x = x + dx;
    % stop condition
    if norm(dx) < threshold
        disp("Converged!");
        break;
    end
end
%% Part 2 ouput x
disp("--- Final Result for Case IV (h1, h4) ---");
fprintf('Total Iterations: %d\n', iter);
disp("Solution [xD, xE, yD, yE] = ");
disp(x);

%% Part 3 ouput Cofactor matrix
A_inv = inv(A);
Q_xx = A_inv(1:4, 1:4);

disp("Cofactor Matrix Q_xx (for [xD, xE, yD, yE]) = 1e-3 *");
disp(Q_xx * 1e3);

%% Part 4 Check constriants (h1, h4)
disp("--- Constraint Verification ---");
xD_final = x(1); xE_final = x(2); yD_final = x(3); yE_final = x(4);

h1_final = -yD_final * xE_final - 4 * yE_final + 4 * yD_final + yE_final * xD_final;
h4_final = 2 * xD_final - xE_final - 4;

fprintf('Final constraint value h1(x) = %.2e\n', h1_final);
fprintf('Final constraint value h4(x) = %.2e\n', h4_final);

tolerance = 1e-12; 
if all(abs([h1_final, h4_final]) < tolerance)
    disp('Conclusion: Constraints h1 and h4 are satisfied.');
else
    disp('Conclusion: The constraints are NOT satisfied.');
end