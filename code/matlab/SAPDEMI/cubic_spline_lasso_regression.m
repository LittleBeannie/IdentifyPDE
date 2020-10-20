function [beta_lasso, y_lasso, X_lasso, fitinfo] = cubic_spline_lasso_regression(uNoise, xMax, xNum, tMax, tNum, lambda, interaction, smoothing, smoothing_alpha)
%% input arguments:
%%% uNoise:  the observed values
%%% xMax: the maximum of the spatial variables
%%% xNum: spatial resolution
%%% tMax: the maximum of the temporal variable
%%% tNum: temporal resolution
%%% lambda: the lasso penalty parameter
%%% interaction: 0 OR 1. 0: without interaction; 1: with interaction
%%% smoothing: 0 OR 1. 0: without smoothing; 1: with smoothing
%%% smoothing_alpha: a postive number between 0 and 1, which trades off the smoothing and fitness

%% output returns:
%%% beta_lasso: the estimated parameter from the LASSO model
%%% y_lasso: the response vector in the LASSO model
%%% X_lasso: the design matrix in the LASSO model
%%% fitinfo:  the fit information in the LASSO model


length_T = size(uNoise,1);
length_X = size(uNoise,2);
dx    = xMax/xNum;
dt    = tMax/tNum;

spline_A = zeros(length_X-2, length_X);   %% to estimate f_hat with fixed t
for i = 1:(length_X-2)
    spline_A(i,i) = 1/dx;
    spline_A(i,i+1) = -2/dx;
    spline_A(i,i+2) = 1/dx;
end

spline_AA = zeros(length_T-2, length_T);  %% to estimate f_hat with fixed x
for i = 1:(length_T-2)
    spline_AA(i,i) = 1/dt;
    spline_AA(i,i+1) = -2/dt;
    spline_AA(i,i+2) = 1/dt;
end

spline_B = zeros(length_X, length_X);     %% to estimate derivative w.r.t x1
spline_B(1,1) = -3/(dx^2);
spline_B(1,2) =  3/(dx^2);
spline_B(length_X,length_X-1) =  -3/(dx^2);
spline_B(length_X,length_X)   =   3/(dx^2);
for i = 2:(length_X-1)
    spline_B(i,i-1) = -3/(dx^2);
    spline_B(i,i) = 0;
    spline_B(i,i+1) = 3/(dx^2);
end


spline_BB = zeros(length_T, length_T);    %% to estimate derivative w.r.t t1
spline_BB(1,1) = -3/(dt^2);
spline_BB(1,2) =  3/(dt^2);
spline_BB(length_T,length_T-1) =  -3/(dt^2);
spline_BB(length_T,length_T)   =   3/(dt^2);
for i = 2:(length_T-1)
    spline_BB(i,i-1) = -3/(dt^2);
    spline_BB(i,i) = 0;
    spline_BB(i,i+1) = 3/(dt^2);
end

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

spline_M_x2 = zeros(length_X-2, length_X-2);
spline_M_x2(1,1) = 2*dx/3;
spline_M_x2(1,2) = dx/6;
spline_M_x2(length_X-2,length_X-2) = 2*dx/3;
spline_M_x2(length_X-2,length_X-3) = dx/6;
for i = 2:length_X-3
    spline_M_x2(i,i) = 2*dx/3;
    spline_M_x2(i,i+1) = dx/6;
    spline_M_x2(i,i-1) = dx/6;
end

spline_M_t2 = zeros(length_T-2, length_T-2);
spline_M_t2(1,1) = 2*dt/3;
spline_M_t2(1,2) = dt/6;
spline_M_t2(length_T-2,length_T-2) = 2*dt/3;
spline_M_t2(length_T-2,length_T-3) = dt/6;
for i = 2:length_T-3
    spline_M_t2(i,i) = 2*dt/3;
    spline_M_t2(i,i+1) = dt/6;
    spline_M_t2(i,i-1) = dt/6;
end


if smoothing == 0
    D_t1 = zeros(length_T, length_T);
    L_t1 = eye(length_T, length_T);
    D_t1(1,1) = spline_Q_t1(1,1);
    for i=1:(length_T-1)
        L_t1(i+1, i)   = spline_Q_t1(i,i+1)/D_t1(i,i);
        D_t1(i+1, i+1) = spline_Q_t1(i+1,i+1) - D_t1(i,i)*L_t1(i+1, i)^2;
    end
    Dinv_t1 = inv(D_t1);
    
    D_x1 = zeros(length_X, length_X);
    L_x1 = eye(length_X, length_X);
    D_x1(1,1) = spline_Q_x1(1,1);
    for i=1:(length_X-1)
        L_x1(i+1, i)   = spline_Q_x1(i,i+1)/D_x1(i,i);
        D_x1(i+1, i+1) = spline_Q_x1(i+1,i+1) - D_x1(i,i)*L_x1(i+1, i)^2;
    end
    Dinv_x1 = inv(D_x1);
    
    D_x2 = zeros(length_X-2, length_X-2);
    L_x2 = eye(length_X-2, length_X-2);
    D_x2(1,1) = 2*dx/3;
    for i=1:(length_X-3)
        L_x2(i+1, i) = dx/6/D_x2(i,i);
        D_x2(i+1, i+1) = 2*dx/3 - dx^2/36/D_x2(i,i);
    end
    Dinv_x2 = inv(D_x2);
    
    %% solve the derivative w.r.t t
    u_derivative_t = [];
    for x_fix = 1: length_X
        uNoise_fix_x = uNoise(:,x_fix);
        spline_y = 3*diff(uNoise_fix_x(1:end-1))/dt/dt + 3*diff(uNoise_fix_x(2:end))/dt/dt ;
        spline_y = [3*(uNoise_fix_x(2)-uNoise_fix_x(1))/dt/dt; spline_y; 3*(uNoise_fix_x(end)-uNoise_fix_x(end-1))/dt/dt];
        spline_theta = (L_t1'\Dinv_t1)*(L_t1\spline_y);
        %spline_theta = spline_Q_t1\spline_y;
        u_derivative_t = [u_derivative_t,spline_theta];
    end
    
    %% solve the first derivative w.r.t. x1
    u_derivative_x1 = [];
    for t_fix = 1: length_T
        uNoise_fix_t = uNoise(t_fix,:)';
        spline_y = 3*diff(uNoise_fix_t(1:end-1))/dx/dx + 3*diff(uNoise_fix_t(2:end))/dx/dx ;
        spline_y = [3*(uNoise_fix_t(2)-uNoise_fix_t(1))/dx/dx; spline_y; 3*(uNoise_fix_t(end)-uNoise_fix_t(end-1))/dx/dx];
        spline_theta = (L_x1'\Dinv_x1)*(L_x1\spline_y);
        %spline_theta = spline_Q_x1\spline_y;
        u_derivative_x1 = [u_derivative_x1, spline_theta];
    end
    
    %% solve the first derivative w.r.t. x2
    u_derivative_x2 = [];
    for t_fix = 1: length_T
        uNoise_fix_t = uNoise(t_fix,:);
        spline_y = ( diff(uNoise_fix_t(2:end))/dx - diff(uNoise_fix_t(1:(end-1)))/dx)';
        spline_sigma = (L_x2'\Dinv_x2)*(L_x2\spline_y);
        %spline_sigma = spline_M_x2\spline_y;
        spline_sigma_extend = [0; spline_sigma; 0];
        u_derivative_x2 =[u_derivative_x2, spline_sigma_extend];
    end
    
    
    %% lasso regression (prior in spatio)
    y_lasso = u_derivative_t'; y_lasso = y_lasso(:);
    uNoise_trans = uNoise';
    if interaction == 0
        X_lasso = [ones(size(u_derivative_x2,1)*size(u_derivative_x2,2),1), uNoise_trans(:),  u_derivative_x1(:) ,     u_derivative_x2(:)    ];
        [beta_lasso, fitinfo] = lasso(X_lasso, y_lasso, 'Lambda', lambda);
    elseif interaction == 1
        X_lasso = [ones(size(u_derivative_x2,1)*size(u_derivative_x2,2),1), uNoise_trans(:),  u_derivative_x1(:) ,     u_derivative_x2(:) ,  uNoise_trans(:).^2,  u_derivative_x1(:).^2,  u_derivative_x2(:).^2, uNoise_trans(:) .* u_derivative_x1(:),  uNoise_trans(:) .* u_derivative_x2(:),   u_derivative_x1(:).* u_derivative_x2(:) ];
        [beta_lasso, fitinfo] = lasso(X_lasso, y_lasso, 'Lambda', lambda);
    end
    
    
elseif smoothing == 1
     %% solve the derivative w.r.t t
     u_derivative_t = [];
     for x_fix = 1: length_X
         uNoise_fix_x = uNoise(:,x_fix);
         spline_f = ( smoothing_alpha * eye(length_T) + (1-smoothing_alpha) * (spline_AA') * spline_M_t2 * spline_AA)\(smoothing_alpha * uNoise_fix_x);
         spline_theta = spline_Q_t1 \(spline_BB * spline_f);
         u_derivative_t = [u_derivative_t, spline_theta];
     end
     %% solve the first derivative w.r.t. x0, x1, x2
     u_derivative_x0 = [];
     u_derivative_x1 = [];
     u_derivative_x2 = [];
     for t_fix = 1: length_T
        uNoise_fix_t = uNoise(t_fix,:)';
        spline_f = (smoothing_alpha*eye(length_X) + (1-smoothing_alpha)* (spline_A') * spline_M_x2 * spline_A)\(smoothing_alpha*uNoise_fix_t);
        spline_theta = spline_Q_x1 \(spline_B * spline_f);
        spline_sigma = spline_M_x2 \(spline_A * spline_f);
        
        u_derivative_x0 = [u_derivative_x0, spline_f];
        u_derivative_x1 = [u_derivative_x1, spline_theta];
        u_derivative_x2 = [u_derivative_x2, [0; spline_sigma; 0]];
    end
     
     %% lasso regression (prior in spatio)
     y_lasso = u_derivative_t'; y_lasso = y_lasso(:);
     if interaction == 0
         X_lasso = [ones(size(u_derivative_x2,1)*size(u_derivative_x2,2),1), u_derivative_x0(:),  u_derivative_x1(:) ,     u_derivative_x2(:)    ];
         [beta_lasso, fitinfo] = lasso(X_lasso, y_lasso, 'Lambda', lambda);
     elseif interaction == 1
         X_lasso = [ones(size(u_derivative_x2,1)*size(u_derivative_x2,2),1), u_derivative_x0(:),  u_derivative_x1(:) ,     u_derivative_x2(:) ,  u_derivative_x0(:).^2,  u_derivative_x1(:).^2,  u_derivative_x2(:).^2, u_derivative_x0(:) .* u_derivative_x1(:),  u_derivative_x0(:) .* u_derivative_x2(:),   u_derivative_x1(:).* u_derivative_x2(:) ];
         [beta_lasso, fitinfo] = lasso(X_lasso, y_lasso, 'Lambda', lambda);
     end
end
    











   

%save('cubic_spline_rideg_regression_Mat')
end