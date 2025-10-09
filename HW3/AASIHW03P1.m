% AASI HW3
% Po-Tao, Lin | B11501037 |
% Date: 2025/10/02

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
rand_err = sigma .* randn(nEdges,1);

meas_dh = true_dh + rand_err;

%% Part 4 compute
H1 = 17.85;
H_diff = zeros(16,1);
H_diff([1,6,7]) = [-H1,H1,-H1];
f = -meas_dh+H_diff;
B = [0,0,0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,-1,0,0,0,0,0;
    0,1,0,0,0,-1,0,0,0,0,0,0;
    0,-1,0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,1,0,-1,0;
    0,0,0,0,0,0,0,0,-1,0,0,0;
    0,0,0,0,1,0,0,0, 0,0,0,0;
    0,0,1,0,-1,0,0,0,0,0,0,0;
    0,0,-1,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,1,0,0,-1,0,0;
    0,0,0,0,0,0,0,1,0,-1,0,0;
    0,0,0,0,0,0,0,-1,0,0,0,1;
    0,0,0,0,0,0,0,0,0,0,1,-1;
    1,0,0,0,0,0,0,0,0,0,0,-1;
    -1,0,0,1,0,0,0,0,0,0,0,0;
    0,0,0,-1,1,0,0,0,0,0,0,0
];
%% Part 4-1 equal weight
W = diag(ones(1,16));           %w, weight matrix
x = inv((B')*W*B)*(-B')*W*f;    %x, elevation of H_2 to H_13
v = B*x+f;
%% Part 4-2 unequal weight

sigma_0 = 0.007.* sqrt(dist_km); %7mm/sqrt(km)
sigma_i = sigma_0;
W = diag(1./sigma_i.^2);        %w, weight matrix
x = inv((B')*W*B)*(-B')*W*f;    %x, elevation of H_2 to H_13
v = B*x+f;

%% Part 5 different sigma 0
% 5-1 ->14mm/sqrt(km)
sigma_0 = 0.014.* sqrt(dist_km); %14mm/sqrt(km)
sigma_i = sigma_0;
W = diag(1./sigma_i.^2);        %w, weight matrix
x1 = inv((B')*W*B)*(-B')*W*f;    %x, elevation of H_2 to H_13
v1 = B*x1+f;
% 5-2->3.5mm/sqrt(km)
sigma_0 = 0.0035.* sqrt(dist_km); %3.5mm/sqrt(km)
sigma_i = sigma_0;
W = diag(1./sigma_i.^2);        %w, weight matrix
x2 = inv((B')*W*B)*(-B')*W*f;    %x, elevation of H_2 to H_13
v2 = B*x2+f;
%% Part 6 different weight scheme
sigma_0 = 0.007.* sqrt(dist_km); %7mm/sqrt(km)
sigma_i = sigma_0;
% 6-1->wi = 1=1/sigma_i
W = diag(1./sigma_i);
x3 = inv((B')*W*B)*(-B')*W*f;    %x, elevation of H_2 to H_13
v3 = B*x3+f;
% 6-2->wi = exp(-100sigma_i)
W = diag(exp(-100*sigma_i));
x4 = inv((B')*W*B)*(-B')*W*f;    %x, elevation of H_2 to H_13
v4 = B*x4+f;