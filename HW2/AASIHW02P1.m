% AASI HW2
% Po-Tao, Lin | B11501037 |
% Date: 2025/09/25

clear; close all; clc;


%% Part 1: generate network 
rng(3);  % fix random seed (convinience to regenerate)

nPts = 13;
coords = 400*rand(nPts,2);      % XY-location [m]
elev   = 20 + 5*randn(nPts,1);  % elevation [m]

% network (13 points, 4 loops)
edges = [1 8; 8 7; 7 3; 3 12; 12 10; 10 1;  % loop 1: [1 8 7 3 12 10]
        1 6; 6 4; 4 11; 11 8;               % loop 2: [1 6 4 11 8]
        11 9; 9 13; 13 12;                  % loop 3
        13 2; 2 5; 5 6;];                   % loop 4            

nEdges = size(edges,1);

%% Part 2 true level difference
true_dh = elev(edges(:,2)) - elev(edges(:,1));

%% Part 3 add random error (7mm/sqrt(km))
dist_km = sqrt(sum((coords(edges(:,2),:) - coords(edges(:,1),:)).^2,2))/1000; % distance [km]
sigma = 0.007.* sqrt(dist_km); %7mm/sqrt(km)


%% ---------matrix form
%% Part 4 for 2-4
f = [0.007;0.002;0;-0.004];
A = [1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0;
    -1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0;
    0,-1,-1,-1,0,0,0,0,0,1,1,1,1,0,0,0;
    0,0,0,0,-1,-1,-1,0,0,0,0,0,-1,1,1,1];
sigma_0 = 0.007.* sqrt(dist_km); %14mm/sqrt(km)
sigma_i = sigma_0;
W = diag(1./sigma_i.^2);
Q = inv(W);
v = Q*A'*((A*Q*A')\f);
%% Part 5 for 2-6
% 2-6-1->14mm/sqrt(km)
sigma_0 = 0.014.* sqrt(dist_km); %14mm/sqrt(km)
sigma_i = sigma_0;
W = diag(1./sigma_i.^2);
Q = inv(W);
v1 = Q*A'*((A*Q*A')\f);
% 2-6-2->3.5mm/sqrt(km)
sigma_0 = 0.0035.* sqrt(dist_km); %3.5mm/sqrt(km)
sigma_i = sigma_0;
W = diag(1./sigma_i.^2);
Q = inv(W);
v2 = Q*A'*((A*Q*A')\f);
% 2-6-3->7mm/km
sigma_0 = 0.007.* dist_km; %7mm/km
sigma_i = sigma_0;
W = diag(1./sigma_i.^2);
Q = inv(W);
v3 = Q*A'*((A*Q*A')\f);
%% Part 5 for 2-7
sigma_0 = 0.007.* sqrt(dist_km); %7mm/sqrt(km)
sigma_i = sigma_0;
% 2-7-1->wi = 1
W = diag(ones(size(sigma_i)));
Q = inv(W);
v4 = Q*A'*((A*Q*A')\f);
% 2-7-2->wi = 1=1/sigma_i
W = diag(1./sigma_i);
Q = inv(W);
v5 = Q*A'*((A*Q*A')\f);
% 2-7-3->wi = exp(-100sigma_i)
W = diag(exp(-100*sigma_i));
Q = inv(W);
v6 = Q*A'*((A*Q*A')\f);