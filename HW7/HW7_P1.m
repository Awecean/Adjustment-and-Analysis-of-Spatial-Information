% AASI HW07_P1
% B11501037 Tao

%% Part 0 clear
clear; clc; close all;
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
x = [8; 12; 4; 8];   % same as your initial X_U_a

%---------------------------
% Iteration parameters
%---------------------------
W = eye(4);
threshold = 1e-10;
%% Part 2 iteration
for iter = 1:20

    xD = x(1);  xE = x(2);
    yD = x(3);  yE = x(4);

    %---------------------------
    % Nonlinear equations F(x)
    %---------------------------
    F = [
        (xE-0)^2 + (yE-4)^2 - s1^2;
        (xD-0)^2 + (yD-4)^2 - s2^2;
        (xD-8)^2 + (yD-0)^2 - s3^2;
        (xE-8)^2 + (yE-0)^2 - s4^2
    ];

    %---------------------------
    % Jacobian matrix B = dF/dx
    %---------------------------
    B = [0,2*xE, 0, 2*yE-8;
        2*xD,0,2*yD-8, 0;
        2*xD-16,0, 2*yD, 0;
        0,2*xE-16, 0, 2*yE];

    %---------------------------
    % Solve for correction: dx
    %---------------------------
    dx = -inv(B' * W * B)*(B' * W *F);

    x = x+dx;   % update

    % stop condition
    if norm(dx) < threshold
        break;
    end
end
%% Part 3 Ouput
disp("Solution [xD, xE, yD, yE] = ");
disp(x);
disp("Cofactor Matrix = 1e-3*")
disp(inv(B'*W*B)*1e3);