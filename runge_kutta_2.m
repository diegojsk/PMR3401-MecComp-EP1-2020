%% PMR 3401 - Mecânica Computacional para Mecatrônica - 2020
%  Exercício Programa 1
%  Implementação do Método de Runge-Kutta de 2 ordem
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 10340624
%
%% M�todo de Rungwe-Kutta de 2 ordem
function [X, Y] = runge_kutta_2(f,h,yi,xi,xf)
% f : fun��o do tipo y' = f(x, y)
% h : passo 
% yi, xi : condi��es iniciais
% xf : condi��o de parada
i = 1;
y(:,i) = yi;
x(:,i) = xi;
steps = (xf-xi)/h;

% La�o Principal 
for i=1:steps
    k1 = f(x(:,i),y(:,i));
    k2 = f(x(:,i) + h*(1/2), y(:,i) +h*k1*(1/2) );
    y(:,i+1) = y(:,i) + h*k2;
    x(:,i+1) = x(:,i) + h;
end
 Y = y;
 X = x;
end