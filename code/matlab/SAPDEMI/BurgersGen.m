%%%%%%%%%% Generate Bateman-Burgers Example Data %%%%%%%%%%%%%%%%
%%%%%%%%%% Roy  Jan.2.2020 %%%%%%%%
%%%% Periodic boundary condition
%%%% Lax-Wendroff Scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = BurgersGen(f,dx,dt,xMax,tMax,nu)
M = floor(xMax/dx);
N = floor(tMax/dt);

U = zeros(N+1,M+1);
U(1,:) = f(0:dx:dx*M);
R = dt/dx;
r = nu*dt/dx^2;
for n = 2:N+1
    U(n,2:end-1) = U(n-1,2:end-1)-...
        R/2*U(n-1,2:end-1).*(U(n-1,3:end)-U(n-1,1:end-2))+...
        (R^2/2*U(n-1,2:end-1).^2+r).*(U(n-1,3:end)-2*U(n-1,2:end-1)+U(n-1,1:end-2));
end
end