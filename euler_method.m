%% PMR 3401 - Mec�nica Computacional para Mecatr�nica - 2020
%  Exerc�cio Programa 1
%  Implementa��o do M�todo de Euler
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 10340624
%
%% Euler's Method

function [X, Y] = euler_method(f,h,yi,xi,xf)
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
 Y = y;
 X = x;
end