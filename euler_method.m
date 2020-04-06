%% PMR 3401 - Mecï¿½nica Computacional para Mecatrï¿½nica - 2020
%  Exercï¿½cio Programa 1
%  Implementaï¿½ï¿½o do Mï¿½todo de Euler
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 10340624
%
%% Método de Euler
function [X, Y] = euler_method(f,h,yi,xi,xf)
% f : função do tipo y' = f(x, y)
% h : passo 
% yi, xi : condições iniciais
% xf : condição de parada

i = 1;
y(:,i) = yi;
x(:,i) = xi;
steps = (xf-xi)/h;

% Laço Principal do Método de Euler
for i=1:steps
    y(:,i+1) = y(:,i) + h*f(x(:,i),y(:,i));
    x(:,i+1) = x(:,i) + h;
end
Y = y;
X = x;
end