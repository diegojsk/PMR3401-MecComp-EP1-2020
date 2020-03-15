%% PMR 3401 - Mec�nica Computacional para Mecatr�nica - 2020
%  Exerc�cio Programa 1
%  Implementa��o do M�todo de Runge-Kutta de 2 ordem
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 
%
%% Runge-Kutta 2 order Method

function [X, Y] = runge_kutta_2(f,h,yi,xi,xf)
%% Algorithm
i = 1;
y(:,i) = yi;
x(:,i) = xi;
steps = (xf-xi)/h;
for i=1:steps
    k1 = f(x(:,i),y(:,i));
    k2 = f(x(:,i) + h*(1/2), y(:,i) +h*k1*(1/2) );
    y(:,i+1) = y(:,i) + h*k2;
    x(:,i+1) = x(:,i) + h;
end
 Y = y;  
 X = x;
end