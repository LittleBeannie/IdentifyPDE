function [beta_ridge, y_ridge, X_ridge] = cubic_spline_ridge_regression(uNoise, xMax, xNum, tMax, tNum, lambda, interaction)
%% input arguments:
%%% uNoise:  the observed values
%%% xMax: the maximum of the spatial variables
%%% xNum: spatial resolution
%%% tMax: the maximum of the temporal variable
%%% tNum: temporal resolution
%%% lambda: the ridged regression penalty parameter
%%% interaction: 0 OR 1. 0: without interaction; 1: with interaction

%% output returns:
%%% beta_lasso: the estimated parameter from the LASSO model
%%% y_lasso: the response vector in the LASSO model
%%% X_lasso: the design matrix in the LASSO model



length_T = size(uNoise,1);
length_X = size(uNoise,2);
dx    = xMax/xNum;
dt    = tMax/tNum;


spline_Q_t1 = zeros(length_T, length_T);
spline_Q_t1(1,1) = 2/dt;
spline_Q_t1(1,2) = 1/dt;
spline_Q_t1(length_T,length_T-1) = 1/dt;
spline_Q_t1(length_T,length_T) = 2/dt;
for i = 2:(length_T-1)
    spline_Q_t1(i,i-1) = 1/dt;
    spline_Q_t1(i,i) = 4/dt;
    spline_Q_t1(i,i+1) = 1/dt;
end

spline_Q_x1 = zeros(length_X, length_X);
spline_Q_x1(1,1) = 2/dx;
spline_Q_x1(1,2) = 1/dx;
spline_Q_x1(length_X,length_X-1) = 1/dx;
spline_Q_x1(length_X,length_X) = 2/dx;
for i = 2:(length_X-1)
    spline_Q_x1(i,i-1) = 1/dx;
    spline_Q_x1(i,i) = 4/dx;
    spline_Q_x1(i,i+1) = 1/dx;
end


% D_t1 = zeros(length_T, length_T);
% Dinv_t1 = zeros(length_T, length_T);
% L_t1 = eye(length_T, length_T);
% D_t1(1,1) = spline_Q_t1(1,1);
% Dinv_t1(1,1) = 1/D_t1(1,1);
% for i=1:(length_T-1)
%     L_t1(i+1, i)   = spline_Q_t1(i,i+1)/D_t1(i,i);
%     D_t1(i+1, i+1) = spline_Q_t1(i+1,i+1) - D_t1(i,i)*L_t1(i+1, i)^2;
%     Dinv_t1(i+1, i+1) = 1/D_t1(i+1, i+1);
% end
% L_t1_inv = -L_t1 + 2* eye(size(L_t1,2));
% for i = 2:size(L_t1, 1)
%     for j = 1: (i -1)
%         L_t1_inv(i,i-j) = -L_t1(i, i-1) * L_t1_inv(i-1, i-j);
%     end
% end

D_t1_diag = zeros(length_T,1);
D_t1_diag(1) = spline_Q_t1(1,1);
L_t1_offdiag = zeros(length_T - 1, 1);
for i=1:(length_T-1)
    L_t1_offdiag(i) = spline_Q_t1(i,i+1)/D_t1_diag(i);
    D_t1_diag(i+1) = spline_Q_t1(i+1,i+1) - D_t1_diag(i) * L_t1_offdiag(i)^2;
end
    


% 
% D_x1 = zeros(length_X, length_X);
% Dinv_x1 = zeros(length_X, length_X);
% L_x1 = eye(length_X, length_X);
% D_x1(1,1) = spline_Q_x1(1,1);
% Dinv_x1(1,1) = 1/D_x1(1,1);
% for i=1:(length_X-1)
%     L_x1(i+1, i)   = spline_Q_x1(i,i+1)/D_x1(i,i);
%     D_x1(i+1, i+1) = spline_Q_x1(i+1,i+1) - D_x1(i,i)*L_x1(i+1, i)^2;
%     Dinv_x1(i+1, i+1) = 1/D_x1(i+1, i+1);
% end
% L_x1_inv = -L_x1 + 2* eye(size(L_x1,2));
% for i = 2:size(L_x1, 1)
%     for j = 1: (i -1)
%         L_x1_inv(i,i-j) = -L_x1(i, i-1) * L_x1_inv(i-1, i-j);
%     end
% end
D_x1_diag = zeros(length_X,1);
D_x1_diag(1) = spline_Q_x1(1,1);
L_x1_offdiag = zeros(length_X - 1, 1);
for i=1:(length_X-1)
    L_x1_offdiag(i) = spline_Q_x1(i,i+1)/D_x1_diag(i);
    D_x1_diag(i+1) = spline_Q_x1(i+1,i+1) - D_x1_diag(i) * L_x1_offdiag(i)^2;
end



% D_x2 = zeros(length_X-2, length_X-2);
% Dinv_x2 = zeros(length_X-2, length_X-2);
% L_x2 = eye(length_X-2, length_X-2);
% D_x2(1,1) = 2*dx/3;
% Dinv_x2(1,1) = 1/D_x2(1,1);
% for i=1:(length_X-3)
%     L_x2(i+1, i) = dx/6/D_x2(i,i);
%     D_x2(i+1, i+1) = 2*dx/3 - dx^2/36/D_x2(i,i);
%     Dinv_x2(i+1, i+1) = 1/D_x1(i+1, i+1);
% end
% L_x2_inv = -L_x2 + 2* eye(size(L_x2,2));
% for i = 2:size(L_x2, 1)
%     for j = 1: (i -1)
%         L_x2_inv(i,i-j) = -L_x2(i, i-1) * L_x2_inv(i-1, i-j);
%     end
% end
D_x2_diag = zeros(length_X-2,1);
D_x2_diag(1) = 2*dx/3;
L_x2_offdiag = zeros(length_X - 3, 1);
for i=1:(length_X-3)
    L_x2_offdiag(i) = dx/6/D_x2_diag(i);
    D_x2_diag(i+1) = 2*dx/3 - dx^2/36/D_x2_diag(i);
end



%% solve the derivative w.r.t t
u_derivative_t = [];
for x_fix = 1: length_X
    uNoise_fix_x = uNoise(:,x_fix);
    spline_y = 3*diff(uNoise_fix_x(1:end-1))/dt/dt + 3*diff(uNoise_fix_x(2:end))/dt/dt ;
    spline_y = [3*(uNoise_fix_x(2)-uNoise_fix_x(1))/dt/dt; spline_y; 3*(uNoise_fix_x(end)-uNoise_fix_x(end-1))/dt/dt];
    
    spline_yy = zeros(length_T, 1);  %% product of L^{-1} * spline_y
    spline_yy(1) = spline_y(1);
    for i = 2: length_T
        spline_yy(i) = spline_y(i) - L_t1_offdiag(i-1)*spline_yy(i-1);
    end
    
    spline_yyy = spline_yy ./ D_t1_diag;  %% product of D^{-1} * L^{-1} * spline_y
    
    spline_yyyy = zeros(length_T, 1);     %% product of L'^{-1}* D^{-1} * L^{-1} * spline_y
    spline_yyyy(length_T) = spline_yyy(length_T);
    ii = length_T-1;
    while ii>0
        spline_yyyy(ii) = spline_yyy(ii) - L_t1_offdiag(ii) * spline_yyyy(ii+1);
        ii = ii-1;
    end
    %spline_theta = L_t1_inv' * Dinv_t1 * L_t1_inv * spline_y;
    u_derivative_t = [u_derivative_t,  spline_yyyy];
end
    
%% solve the first derivative w.r.t. x1
u_derivative_x1 = [];
for t_fix = 1: length_T
    uNoise_fix_t = uNoise(t_fix,:)';
    spline_y = 3*diff(uNoise_fix_t(1:end-1))/dx/dx + 3*diff(uNoise_fix_t(2:end))/dx/dx ;
    spline_y = [3*(uNoise_fix_t(2)-uNoise_fix_t(1))/dx/dx; spline_y; 3*(uNoise_fix_t(end)-uNoise_fix_t(end-1))/dx/dx];
    
    spline_yy = zeros(length_X, 1);  %% product of L^{-1} * spline_y
    spline_yy(1) = spline_y(1);
    for i = 2: length_X
        spline_yy(i) = spline_y(i) - L_x1_offdiag(i-1)*spline_yy(i-1);
    end
    
    spline_yyy = spline_yy ./ D_x1_diag;  %% product of D^{-1} * L^{-1} * spline_y
    
    spline_yyyy = zeros(length_X, 1);     %% product of L'^{-1}* D^{-1} * L^{-1} * spline_y
    spline_yyyy(length_X) = spline_yyy(length_X);
    ii = length_X-1;
    while ii>0
        spline_yyyy(ii) = spline_yyy(ii) - L_x1_offdiag(ii) * spline_yyyy(ii+1);
        ii = ii-1;
    end
    
    %spline_theta = L_x1_inv' * Dinv_x1 * L_x1_inv * spline_y;
    u_derivative_x1 = [u_derivative_x1, spline_yyyy];
end
    
%% solve the first derivative w.r.t. x2
u_derivative_x2 = [];
for t_fix = 1: length_T
    uNoise_fix_t = uNoise(t_fix,:);
    spline_y = ( diff(uNoise_fix_t(2:end))/dx - diff(uNoise_fix_t(1:(end-1)))/dx)';
    
    spline_yy = zeros(length_X-2, 1);  %% product of L^{-1} * spline_y
    spline_yy(1) = spline_y(1);
    for i = 2: (length_X-2)
        spline_yy(i) = spline_y(i) - L_x2_offdiag(i-1)*spline_yy(i-1);
    end
    
    spline_yyy = spline_yy ./ D_x2_diag;  %% product of D^{-1} * L^{-1} * spline_y
    
    spline_yyyy = zeros(length_X-2, 1);     %% product of L'^{-1}* D^{-1} * L^{-1} * spline_y
    spline_yyyy(length_X-2) = spline_yyy(length_X-2);
    ii = length_X-3;
    while ii>0
        spline_yyyy(ii) = spline_yyy(ii) - L_x2_offdiag(ii) * spline_yyyy(ii+1);
        ii = ii-1;
    end
    %spline_sigma = L_x2_inv' * Dinv_x2 * L_x2_inv * spline_y;
    u_derivative_x2 =[u_derivative_x2, [0; spline_yyyy; 0]];
end
    








%% ridge regression (prior in spatio)
y_ridge = u_derivative_t'; y_ridge = y_ridge(:);
uNoise_trans = uNoise';
if interaction == 0
    X_ridge = [ones(size(u_derivative_x2,1)*size(u_derivative_x2,2),1), uNoise_trans(:),  u_derivative_x1(:) ,     u_derivative_x2(:)    ];
    beta_ridge = (X_ridge' * X_ridge + lambda * eye(size(X_ridge',1)))\(X_ridge' *y_ridge);
    save('info')
elseif interaction == 1
    X_ridge = [ones(size(u_derivative_x2,1)*size(u_derivative_x2,2),1), uNoise_trans(:),  u_derivative_x1(:) ,     u_derivative_x2(:) ,  uNoise_trans(:).^2,  u_derivative_x1(:).^2,  u_derivative_x2(:).^2, uNoise_trans(:) .* u_derivative_x1(:),  uNoise_trans(:) .* u_derivative_x2(:),   u_derivative_x1(:).* u_derivative_x2(:) ];
    beta_ridge = (X_ridge' * X_ridge + lambda * eye(size(X_ridge',1)))\(X_ridge' *y_ridge);
    save('info')
end
   



end