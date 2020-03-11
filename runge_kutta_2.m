function Y = runge_kutta_2(f,h,yi,xi,xf)
%% Runge-Kutta 2 order Method
%
%% Algorithm
i = 1;
y(:,i) = yi;
x(:,i) = xi;
steps = (xf-xi)/h;
for i=1:steps
    k1 = f(x(:,i),y(:,i));
    k2 = f(x(:,i) + h*(1/2), y(i) +h*k1*(1/2) );
    y(:,i+1) = y(:,i) + h*k2;
    x(:,i+1) = x(:,i) + h;
end
 Y =  y;  

end