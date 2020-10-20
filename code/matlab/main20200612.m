%%% ---------------------------------------------------------------%%
%%%                       SIMULATION 1                             %%
%%%           compare the computational complexity                 %%
%%%            of cubic spline & local polynomial                  %%
%%% ---------------------------------------------------------------%%
addpath(genpath('CountFLOP'))
addpath(genpath('SAPDEMI'))

clear all
clc
sigma = 0.01;
xNum  = 20;
tNum  = 200;
xMax  = 1; 
tMax  = 0.1;
[uNoise, uTure, xData, tData, length_X, length_T] = PDE_data_generation(xNum, tNum, xMax, tMax, sigma, 'transport');
lambda = 1e-2;
%%% Calculate the computational complexity of our method
profile on
[beta_ridge, y_ridge, X_ridge] =  cubic_spline_ridge_regression(uNoise, xMax, xNum, tMax, tNum, lambda, 0);
profileStruct = profile('info');
[flopTotal, Details] = FLOPS('cubic_spline_ridge_regression', 'info.mat', profileStruct);
%%% Calculate the computational complexity of Namjoon's method
profile on
[beta_ridge, y_ridge, X_ridge] =  local_poly_ridge_regression(uNoise, xMax, xNum, tMax, tNum, lambda, 0);
profileStruct = profile('info');
[flopTotal, Details] = FLOPS('local_poly_ridge_regression', 'info.mat', profileStruct);



%%% ---------------------------------------------------------------%%
%%%                       SIMULATION 2                             %%
%%%                  rocovery the support set                      %%
%%% ---------------------------------------------------------------%%
%%% ---------------------------------------------------------------%%
%%%                Example 1: transport equation                   %%
%%% ---------------------------------------------------------------%%

clear all
clc
sigma = 0.01;
xNum  = 100;
tNum  = 40;
xMax  = 1; 
tMax  = 0.1;
lambda = linspace(1e-5,10,100);
[uNoise, uTure, xData, tData, length_X, length_T] = PDE_data_generation(xNum, tNum, xMax, tMax, sigma, 'transport');
[beta_lasso, y_lasso, X_lasso, fitinfo] = cubic_spline_lasso_regression(uNoise, xMax, xNum, tMax, tNum, lambda, 1, 1, 0.8);
%%% plot the LASSO estimation results
figure
for i = 1:size(beta_lasso,1)
    if i == 3
        h1 = semilogx(fitinfo.Lambda, beta_lasso(i,:),'LineWidth',2,'Color','r');
    %elseif i == 8
    %     h2 = semilogx(fitinfo.Lambda,beta_lasso(i,:),'LineWidth',2,'Color','b');
    else    
        semilogx(fitinfo.Lambda, beta_lasso(i,:),'-.','Color','k','LineWidth',1.5)
    end
    hold on
end
xlabel('$\lambda$','Interpreter','Latex')
ylabel('Coefficient','Interpreter','Latex')
hleg1 = legend([h1],{'$u_{x}$'},'Interpreter','Latex');
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',30)


%%% ---------------------------------------------------------------%%
%%%           Example 2: inviscid Burgers' equation                %%
%%% ---------------------------------------------------------------%%
clear all
clc
sigma = 0.01;  %% you can change 'sigma' into 0.01, 0.1 or 1 to generate the following plots
xNum  = 50;
tNum  = 50;
xMax  = 1; 
tMax  = 0.1;
lambda = linspace(1e-2,10,100);
[uNoise, uTure, xData, tData, length_X, length_T] = PDE_data_generation(xNum, tNum, xMax, tMax, sigma, 'inviscid Burger');
[beta_lasso, y_lasso, X_lasso, fitinfo] = cubic_spline_lasso_regression(uNoise, xMax, xNum, tMax, tNum, lambda, 1, 1, 0.8);
%%% plot the LASSO estimation results
figure
for i = 1:size(beta_lasso,1)
    if i == 8
        h1 = semilogx(fitinfo.Lambda, beta_lasso(i,:),'LineWidth',2,'Color','r');
    %elseif i == 8
    %     h2 = semilogx(fitinfo.Lambda,beta_lasso(i,:),'LineWidth',2,'Color','b');
    else    
        semilogx(fitinfo.Lambda, beta_lasso(i,:),'-.','Color','k','LineWidth',1.5)
    end
    hold on
end
xlabel('$\lambda$','Interpreter','Latex')
ylabel('Coefficient','Interpreter','Latex')
hleg1 = legend([h1],{'$uu_{x}$'},'Interpreter','Latex');
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',30)

%%% ---------------------------------------------------------------%%
%%%             Example 3: viscous Burgers' equation               %%
%%% ---------------------------------------------------------------%%
clear all
clc
sigma = 0.01; %% you can change 'sigma' into 0.01, 0.1 or 1 to generate the following plots
xNum  = 50;
tNum  = 50;
xMax  = 1; 
tMax  = 0.1;
lambda = linspace(1e-2,10,100);
[uNoise, uTure, xData, tData, length_X, length_T] = PDE_data_generation(xNum, tNum, xMax, tMax, sigma, 'viscous Burger', 0.1);
[beta_lasso, y_lasso, X_lasso, fitinfo] = cubic_spline_lasso_regression(uNoise, xMax, xNum, tMax, tNum, lambda, 1, 1, 0.8);
%%% plot the LASSO estimation results
figure
for i = 1:size(beta_lasso,1)
    if i == 4
        h1 = semilogx(fitinfo.Lambda, beta_lasso(i,:),'LineWidth',2,'Color','r');
    elseif i == 8
        h2 = semilogx(fitinfo.Lambda,beta_lasso(i,:),'LineWidth',2,'Color','b');
    else    
        semilogx(fitinfo.Lambda, beta_lasso(i,:),'-.','Color','k','LineWidth',1.5)
    end
    hold on
end
xlabel('$\lambda$','Interpreter','Latex')
ylabel('Coefficient','Interpreter','Latex')
hleg1 = legend([h1, h2],{'$u_{xx}$', '$uu_{x}$'},'Interpreter','Latex');
set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',30)




%%% ---------------------------------------------------------------%%
%%%                     data visualization                         %%
%%% ---------------------------------------------------------------%%
%%% ---------------------------------------------------------------%%
%%%                Example 1: transport equation                   %%
%%% ---------------------------------------------------------------%%
clear all
clc
sigma = 0.01; %% you can change 'sigma' into 0.01 or 0.1 to generate the following plots
xNum  = 100;
tNum  = 40;
xMax  = 1; 
tMax  = 0.1;
lambda = linspace(1e-5,10,100);
[uNoise, uTrue, xData, tData, length_X, length_T] = PDE_data_generation(xNum, tNum, xMax, tMax, sigma, 'transport');
[beta_lasso, y_lasso, X_lasso, fitinfo] = cubic_spline_lasso_regression(uNoise, xMax, xNum, tMax, tNum, lambda, 1, 1, 0.8);
uDenoise = reshape(X_lasso(:,2), [xNum+1, tNum]); 
uDenoise = uDenoise';
plotPDE(xData, tNum, uTrue, '$u^*(x,\cdot)$')       %%%% true data
plotPDE(xData, tNum, uNoise, '$u(x,\cdot)$')        %%%% noised data
plotPDE(xData, tNum, uDenoise, '$\hat{u}(x,\cdot)$')%%%% noised data


%%% ---------------------------------------------------------------%%
%%%           Example 2: inviscid Burgers' equation                %%
%%% ---------------------------------------------------------------%%
clear all
clc
sigma = 0.1; %% you can change 'sigma' into 0.01 or 0.1 to generate the following plots
xNum  = 50;
tNum  = 50;
xMax  = 1; 
tMax  = 0.1;
lambda = linspace(1e-5,10,100);
[uNoise, uTrue, xData, tData, length_X, length_T] = PDE_data_generation(xNum, tNum, xMax, tMax, sigma, 'inviscid Burger');
[beta_lasso, y_lasso, X_lasso, fitinfo] = cubic_spline_lasso_regression(uNoise, xMax, xNum, tMax, tNum, lambda, 1, 1, 0.8);
uDenoise = reshape(X_lasso(:,2), [xNum+1, tNum]); 
uDenoise = uDenoise';
plotPDE(xData, tNum, uTrue, '$u^*(x,\cdot)$')       %%%% true data
plotPDE(xData, tNum, uNoise, '$u(x,\cdot)$')        %%%% noised data
plotPDE(xData, tNum, uDenoise, '$\hat{u}(x,\cdot)$')%%%% noised data





%%% ---------------------------------------------------------------%%
%%%             Example 2: viscous Burgers' equation               %%
%%% ---------------------------------------------------------------%%
clear all
clc
sigma = 0.1; %% you can change 'sigma' into 0.01 or 0.1 to generate the following plots
xNum  = 50;
tNum  = 50;
xMax  = 1; 
tMax  = 0.1;
lambda = linspace(1e-5,10,100);
[uNoise, uTrue, xData, tData, length_X, length_T] = PDE_data_generation(xNum, tNum, xMax, tMax, sigma, 'viscous Burger', 0.1);
[beta_lasso, y_lasso, X_lasso, fitinfo] = cubic_spline_lasso_regression(uNoise, xMax, xNum, tMax, tNum, lambda, 1, 1, 0.8);
uDenoise = reshape(X_lasso(:,2), [xNum+1, tNum]); 
uDenoise = uDenoise';
plotPDE(xData, tNum, uTrue, '$u^*(x,\cdot)$')       %%%% true data
plotPDE(xData, tNum, uNoise, '$u(x,\cdot)$')        %%%% noised data
plotPDE(xData, tNum, uDenoise, '$\hat{u}(x,\cdot)$')%%%% noised data
