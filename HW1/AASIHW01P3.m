% AASI HW1 Problem 3
% Po-Tao, Lin | B11501037 |
% Date: 2025/09/18

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

%% Part 3 add random error (±7mm/sqrt(km))
dist_km = sqrt(sum((coords(edges(:,2),:) - coords(edges(:,1),:)).^2,2))/1000; % distance [km]
sigma = 0.007.* sqrt(dist_km); %7mm/sqrt(km)
rand_err = sigma .* randn(nEdges,1);

meas_dh = true_dh + rand_err;

%% Part4 true error
true_err = sigma .* randn(nEdges,1); % generate an true error

%% Part5 ouput ASCII
BS_Pt = edges(:,1);
FS_Pt = edges(:,2); 

% sufficient digits (3)
fmt1 = '%.3f'; % 3 digit right '.' , save to [mm]
fmt2 = '%.2e'; % 2 digit science record (1 is left to the '.')
T_out = table( ...
    BS_Pt, ...
    FS_Pt, ...
    arrayfun(@(x) sprintf(fmt1, x), true_dh, 'UniformOutput', false), ...
    arrayfun(@(x) sprintf(fmt1, x), meas_dh, 'UniformOutput', false), ...
    arrayfun(@(x) sprintf(fmt2, x), true_err, 'UniformOutput', false), ...
    'VariableNames', {'BS_Pt','FS_Pt','dh_true','dh_meas','terror'} ...
);
disp(T_out)
writetable(T_out,'level_network_results.txt','Delimiter','\t','QuoteStrings',false)
%% Part6 Configuration
labels = arrayfun(@(i) sprintf('P%d, H=%.2f m', i, elev(i)), ...
                  1:nPts, 'UniformOutput', false);
figure('Position',[100 100 700 700]); 
set(gcf,'color','white')
hold on; grid on; axis equal
plot(coords(:,1), coords(:,2), 'ko','MarkerFaceColor','y')
textshiftx = 4*[2,2,-4,2,3,2,2,0,-4,3,-13,3,3];
textshifty = 4*[2,2,3,-1,0,0,1,-3,5,0,-4,0,2];
text(coords(:,1)+textshiftx(:), coords(:,2)+textshifty(:), labels,'FontSize',10)
text(120,130,'$\mathrm{I}$','Color','r','Interpreter','latex','FontSize',14)
text(200,70,'$\mathrm{II}$','Color','r','Interpreter','latex','FontSize',14)
text(50,230,'$\mathrm{III}$','Color','r','Interpreter','latex','FontSize',14)
text(320,170,'$\mathrm{IV}$','Color','r','Interpreter','latex','FontSize',14)

for i=1:nEdges %edge
    plot(coords(edges(i,:),1), coords(edges(i,:),2), 'b-')
end

xlim([-80,500]); ylim([-80,400]);
title('Level Network (13 points, 4 loops)')
xlabel('X [m]'); ylabel('Y [m]')
exportgraphics(gca,'P3levelnetwork.pdf');
