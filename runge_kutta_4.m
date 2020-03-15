%% PMR 3401 - Mecânica Computacional para Mecatrônica - 2020
%  Exercício Programa 1
%  Implementação do Método de Runge-Kutta de 4 ordem
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 
%
%% Runge-Kutta 4 order Method

function [X, Y] = runge_kutta_4(f,h,yi,xi,xf)
%% Algorithm
i = 1;
y(:,i) = yi;
x(:,i) = xi;
steps = (xf-xi)/h;
for i=1:steps
    k1 = f(x(:,i),y(:,i));
    k2 = f(x(:,i) + h*(1/2), y(:,i) +h*k1*(1/2) );
    k3 = f(x(:,i) + h*(1/2), y(:,i) +h*k2*(1/2) );
    k4 = f(x(:,i) + h, y(:,i) +h*k3 );
    y(:,i+1) = y(:,i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    x(:,i+1) = x(:,i) + h;
end
 Y = y;  
 X = x;
end