% AASI HW0
% Po-Tao, Lin | B11501037 |
% Date: 2025/09/10

clear; close all; clc;
%% Part 1 Data setting
% original data
A = 20:5:80;                    %age
D = [0.18, 0.12, 0.07, 0.13, 0.21, 0.25, 0.22, 0.37, 0.50,...
    0.68, 0.79, 0.90, 0.96];    %disease rate
%---make figure
figure('Position',[100 100 600 400])
set(gcf, 'Color','White');
scatter(A,D,'k','filled','DisplayName','data')
xlabel('Age (A)');
ylabel('Disease Rate (D)');
xlim([15,85]); ylim([0,1])
exportgraphics(gca,'Fig1.pdf','ContentType','vector')
%% Part 2 Model

%% Part 2-1 polynomial pattern
B1 = [A.^2; A; A.^0];
C1 = D/B1;
A_fit = linspace(min(A), max(A), 200);
model_poly = @(C,x) C(1)*x.^2+C(2)*x+C(3);

%---make figure
figure('Position',[100 100 600 400])
set(gcf, 'Color','White');
scatter(A,D,'k','filled','DisplayName','data')
hold on
xlabel('Age (A)');
ylabel('Disease Rate (D)');
plot(A_fit, model_poly(C1,A_fit),'linewidth',2,'DisplayName','M_1');
title('$M_1$','Interpreter','latex');
exportgraphics(gca,'Fig2.pdf','ContentType','vector');

%% Part 2-2 slowing model
model_exp = @(p, x) p(3) ./ (1 + p(2) * exp(p(1) * x)); % assumed form
p0 = [-0.1, 1, 1];   % first suggestion
% [lsqcurvefit] to fit curve
opts = optimoptions('lsqcurvefit','Display','iter');
lb = [-Inf, 0, 0];      % upperbound (a,b,c)
ub = [Inf, Inf, Inf];   % lowerbound (a,b,c)
p = lsqcurvefit(model_exp, p0, A, D, lb, ub, opts);

A_fit = linspace(min(A), max(A), 200);
D_fit = model_exp(p, A_fit);
%---make figure
figure('Position',[100 100 600 400])
set(gcf, 'Color','White');
scatter(A,D,'k','filled','DisplayName','data')
hold on
xlabel('Age (A)');
ylabel('Disease Rate (D)');
plot(A_fit, D_fit,'b--','linewidth',2,'DisplayName','M_2')
title('$M_2$','Interpreter','latex');
exportgraphics(gca,'Fig3.pdf','ContentType','vector');

%% Part 3 determine the goodness of fit
%--R^2
function RS = RScompute(D_model, D_origin)
    SS_res = sum((D_origin-D_model).^2);
    SS_tot = sum((D_origin-mean(D_origin)).^2);
    RS = 1-SS_res/SS_tot;
end
%--RMSE
function RMSE = RMSEcompute(D_model, D_origin)
    RMSE = sqrt(mean((D_model-D_origin).^2));
end
%
RS_poly = RScompute(model_poly(p,A),D);
RS_exp = RScompute(model_exp(p,A),D);

RMSE_poly = RMSEcompute(model_poly(p,A),D);
RMSE_exp = RMSEcompute(model_exp(p,A),D);

