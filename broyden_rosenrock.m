%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sahebeh Dadboud: 1569395
%sheet 4 - exe 5
%Minimize the Rosenbrock function in 2d in broyden 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%broyden 
function [x,iter,flag] = broyden_rosenrock(f,x0,Jfx0,itermax,tolf,tolx)
 
f = @(x)100*(x(2)-x(1)^2)^2 + (x(1)-1)^2
%x0 = [1.1,1.1]; 
%x0 = [0,0];
x0 = [1.1 , -10];
%testfunction
df = @(x)400*x(1)^3 + (2-400*x(2))*x(1) - 2;
Jfx0 = df(x0);
itermax = 100;
tolf = 1e-7;
tolx = 1e-7;
flag = 0;


if isempty(Jfx0)
    Jfx0 = finite_diff(f,x0);
end

Jfxkinv = Jfx0^(-1); % Inverse exact !!!

changex = 1;

x = zeros(length(x0),itermax);
x(:,1) = x0;

fx = f(x(:,1));

iter = 1;

while (iter < itermax) && (changex > tolx) && (max(abs(fx)) > tolf)
    
    xold = x(:,iter); fxold = fx;
    
    d    = -Jfxkinv*fxold;
    
    x(:,iter+1) = xold + d; 
    
    skp1 = d; 
    
    fx   = f(x(:,iter+1));
    
    ykp1 = fx - fxold;
    
    Jfxkinv = Jfxkinv + (skp1-Jfxkinv*ykp1)/(skp1'*Jfxkinv*ykp1)*(skp1'*Jfxkinv);
    
    changex = max(abs(xold-x(:,iter+1))); %/max(abs(xold));
    
    iter = iter + 1;
end

x = x(:,1:iter);

if (iter >= itermax) 
    flag = 1;
end

if (changex <= tolx) 
    flag = 2;
end

if (max(abs(fx)) <= tolf)
    flag = 3;
end


end


function dG = finite_diff(G,x)
   sz = length(x);
   Gx = G(x);
   dG = zeros(sz,sz);
   epsilon = 1e-8;
   for k = 1:sz 
       y = x; y(k) = y(k) + epsilon;
       Gy = G(y);
       for j = 1:sz
          dG(j,k) = (Gy(j) - Gx(j))/epsilon;
       end
   end
end




