%% PMR 3401 - Mec�nica Computacional para Mecatr�nica - 2020
%  Exerc�cio Programa 1
%  Execu��o de testes para a implementa��o principal
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 
%
clear all
clc
%% Teste
%% Testar Execu��o dos M�todos Num�ricos: 

% y' = 2x => y = x^2
t = @(x,y) 2*x;

euler_method(t,0.1,1,1,5);
runge_kutta_2(t,0.1,1,1,5);
runge_kutta_4(t,0.1,1,1,5);

%% Testar M�todos Num�ricos com entrada function

% @veiculo : function em veiculo.m

euler_method(@veiculo, 0.1,[0; 0.4; 0; -0.1], 0, 60);
runge_kutta_2(@veiculo, 0.1,[0; 0.4; 0; -0.1], 0, 60);
runge_kutta_4(@veiculo, 0.1,[0; 0.4; 0; -0.1], 0, 60);

%% Testar Loop em lista e Plotagem de gr�ficos

i = 1;
for h = [0.1, 0.5, 1]
    [X, Y] = euler_method(t,h,1,1,5);
    figure(i)
    plot(X, Y)
    grid
    str = {'y'};
    legend(str);
    i = i + 1;
end
