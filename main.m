%% PMR 3401 - Mec�nica Computacional para Mecatr�nica - 2019
%  Exerc�cio Programa 1
%
%  Autores: Diego Jun Sato Kurashima
%           Felipe Gomes de Melo D'Elia
%
%% Teste
f = @(x,y) x;
euler_method(f,1,0,0,5)
runge_kutta_2(f,1,0,0,5)
runge_kutta_4(f,1,0,0,5)
