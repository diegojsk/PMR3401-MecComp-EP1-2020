%% PMR 3401 - Mec�nica Computacional para Mecatr�nica - 2020
%  Exerc�cio Programa 1
%  Implementa��o do M�todo de Euler
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 10340624
%
%% M�todo de Euler
function [X, Y] = euler_method(f,h,yi,xi,xf)
% f : fun��o do tipo y' = f(x, y)
% h : passo 
% yi, xi : condi��es iniciais
% xf : condi��o de parada

i = 1;
y(:,i) = yi;
x(:,i) = xi;
steps = (xf-xi)/h;

% La�o Principal do M�todo de Euler
for i=1:steps
    y(:,i+1) = y(:,i) + h*f(x(:,i),y(:,i));
    x(:,i+1) = x(:,i) + h;
end
Y = y;
X = x;
end