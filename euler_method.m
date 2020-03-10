function Y = euler_method(f,h,yi,xi,xf)
%% Euler Method
% Return the solution of a 
%
%%
i = 0;
y(i) = yi;
x(i) = xi;
steps = (xf-xi)/h;
for i=1:steps
    y(i+1) = y(i) + h*f(x(i),y(i));
    x(i+1) = x(i) + h;
end
 Y =  y(i);  
end