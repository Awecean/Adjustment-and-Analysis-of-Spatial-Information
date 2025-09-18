% AASI HW1 Problem 1
% Po-Tao, Lin | B11501037 |
% Date: 2025/09/18

clear; close all; clc;
%% Part 1 Data setting function
function [x,y,z,mean_vals,std_vals] = randn3d(Mx,My,Mz, sx,sy,sz,N)
    x = Mx + sx*randn(N,1);
    y = My + sy*randn(N,1);
    z = Mz + sz*randn(N,1);

    %% (b) Compute mean and std
    mean_vals = [mean(x), mean(y), mean(z)];    % mean value
    std_vals  = [std(x), std(y), std(z)];       % standard deviation

    figure('Position',[100,100,600,400]);
    set(gcf,'Color','White')
    if  N==30
        scatter3(x,y,z,'k','filled','DisplayName','Random Point','AlphaVariable',0.5);
    else
        scatter3(x,y,z,0.3*log((20)),'k','filled','DisplayName','Random Point','AlphaVariable',0.5);
    end

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title([sprintf('%d',N), ' Random Points in 3D']);
    grid on; axis tight manual;
    
    
    
    %% (c) Mark 1σ range
    hold on;
    plot3([mean(x)-std(x), mean(x)+std(x)], [mean(y), mean(y)], [mean(z), mean(z)], 'r-', 'LineWidth',2,'DisplayName',['s_x = ',sprintf('%.3f',std(x))]);
    plot3([mean(x), mean(x)], [mean(y)-std(y), mean(y)+std(y)], [mean(z), mean(z)], 'g-', 'LineWidth',2,'DisplayName',['s_y = ',sprintf('%.3f',std(y))]);
    plot3([mean(x), mean(x)], [mean(y), mean(y)], [mean(z)-std(z), mean(z)+std(z)], 'b-', 'LineWidth',2,'DisplayName',['s_z = ',sprintf('%.3f',std(z))]);
    fprintf('%.3f,%.3f,%.3f',std(x),std(y),std(z))
    legend('Location','northwest')
    exportgraphics(gca,['P1point',sprintf('%d',N),'.pdf'])
end
%% Part 3 Main procedure - Interface
MEANset = zeros(5,3);STDset = zeros(5,3);
Nset = [30,3000,30000,300000,3000000];
for i = 1:5
    
    [~,~,~,MEANset(i,:),STDset(i,:)] =  randn3d(200,100,500,0.2,0.3,0.1,Nset(i));
end
%%
figure('Position',[100,100,600,400])
set(gcf,'Color','White')
scatter(log10(Nset), MEANset(:,1), 'r^','filled');
hold on 
scatter(log10(Nset), MEANset(:,2), 'gs','filled');
scatter(log10(Nset), MEANset(:,3), 'bO','filled');
xticks(log10(Nset));                        
xticklabels(compose('log_{10}(%d)', Nset)); 
xlabel('log_{10} (Point quantity)');ylabel('Sample Mean Value')
yline([200,100,500],'--')
legend('$\bar{x}$','$\bar{y}$','$\bar{z}$','location','eastoutside','Interpreter','latex','Fontsize',12)
exportgraphics(gca,'P1relationMEAN.pdf')


figure('Position',[100,100,600,400])
set(gcf,'Color','White')
scatter(log10(Nset), STDset(:,1), 'r^','filled');
hold on 
scatter(log10(Nset), STDset(:,2), 'gs','filled');
scatter(log10(Nset), STDset(:,3), 'bO','filled');
xticks(log10(Nset));                        
xticklabels(compose('log_{10}(%d)', Nset)); 
xlabel('log_{10} (Point quantity)');ylabel('Sample Standard Deviation')
yline([0.2,0.3,0.1],'--')
legend('$s_x$','$s_y$','$s_z$','location','eastoutside','Interpreter','latex','Fontsize',12)
exportgraphics(gca,'P1relationSTD.pdf')
