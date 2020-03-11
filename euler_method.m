function Y = euler_method(f,h,yi,xi,xf)
%% Euler's Method
% Return the solution of Y of Y' = f(X,Y)
%
% input f : function f of Y' = f(X,Y) 
%       h : step
%     
%% Algorithm
i = 1;
y(:,i) = yi;
x(:,i) = xi;
steps = (xf-xi)/h;
for i=1:steps
    y(:,i+1) = y(:,i) + h*f(x(:,i),y(:,i));
    x(:,i+1) = x(:,i) + h;
end
 Y =  y;  
end