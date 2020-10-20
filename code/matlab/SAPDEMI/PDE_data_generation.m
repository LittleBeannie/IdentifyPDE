function [uNoise,uTrue xData, tData, length_X, length_T] = PDE_data_generation(xNum, tNum, xMax, tMax, sigma, PDEtype, nu)

%% input arguments:
%%% xNum: spatial resolution
%%% tNum: temporal resolution
%%% xMax: the maximum of the spatial variables
%%% tMax: the maximum of the temporal variable
%%% sigma: the Gaussian noise level
%%% PDEtype: the type of PDE where the data is generated, which can be
%%%          chosen form `transport` OR `inviscid Burger` OR `viscous Burger`
%%% nu: a paramter when the `PDEtype` is set as `viscous Burger`

%% output returns:
%%% uNoise: the generated data
%%% uTrue: the ture data
%%% xData: the vector contains all spatial variables
%%% tData:  the vector contains all temporal variables
%%% length_X: number of spatial variables
%%% length_T: number of temporal variables


dx    = xMax/xNum;
dt    = tMax/tNum;
xData = 0:dx:xMax ;
tData = 0:dt:tMax;

if (PDEtype == "transport")
    U = zeros(tNum, xNum+1);
    for n = 1: tNum
        for i = 1:(xNum+1)
            U(n,i) = 2*sin(4*xData(i)-8*tData(n));
        end
    end
    uTrue = U;
    uNoise = U + normrnd(0, sigma, size(U));
    length_T = size(uNoise,1);
    length_X = size(uNoise,2);
    
    
elseif (PDEtype == "inviscid Burger")
    [xMesh,tMesh] = meshgrid(xData,tData); % col dim: x; row dim: t
    fineRatioX = 1; 
    fineRatioT = 5000;
    fineDx = dx/fineRatioX;
    fineDt = dt/fineRatioT;
    fineXData = 0:fineDx:xMax;
    fineTData = 0:fineDt:tMax;
    f = @(x) sin(2*pi*x);                   %(inviscid Burger's equation)
    nu = 0;
    U = BurgersGen(f,fineDx,fineDt,xMax,tMax,nu);
    if fineRatioX==1
        U = U(mod(1:length(fineTData)-1,fineRatioT)==1,:);
    else
        U = U(mod(1:length(fineTData)-1,fineRatioT)==1,...
            mod(1:length(fineXData),fineRatioX)==1);
    end
    uTrue = U;
    rng(3)
    uNoise = U + normrnd(0,sigma,size(U));
    length_T = size(uNoise,1);
    length_X = size(uNoise,2);
    
    
    
elseif (PDEtype == "viscous Burger")
    [xMesh,tMesh] = meshgrid(xData,tData); % col dim: x; row dim: t
    fineRatioX = 1; 
    fineRatioT = 5000;
    fineDx = dx/fineRatioX;
    fineDt = dt/fineRatioT;
    fineXData = 0:fineDx:xMax;
    fineTData = 0:fineDt:tMax;
    f = @(x) sin(4*pi*x).^2+sin(2*pi*x).^3; %(viscous Burger's equation)
    % nu = 0.1;
    U = BurgersGen(f,fineDx,fineDt,xMax,tMax,nu);
    if fineRatioX==1
        U = U(mod(1:length(fineTData)-1,fineRatioT)==1,:);
    else
        U = U(mod(1:length(fineTData)-1,fineRatioT)==1,...
            mod(1:length(fineXData),fineRatioX)==1);
    end
    uTrue = U;
    rng(3)
    uNoise = U + normrnd(0,sigma,size(U));
    length_T = size(uNoise,1);
    length_X = size(uNoise,2);
end
   







end