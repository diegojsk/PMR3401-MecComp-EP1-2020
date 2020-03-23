%% PMR 3401 - Mecânica Computacional para Mecatrônica - 2020
%  Exercício Programa 1
%  Script para analisar solucoes referentes a Questao 2a
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 10340624
%
clear all
clc
%% Questao 2a

L1 = 2;
L2 = 2.5;
L2eixo = 1.8;
m1 = 450;
m2 = 1000;
uIz = 2.7;
R = 0.3;
g = 9.81;
veld = 80/3.6;
F1 = -0.5*m1*g;
F2 = -0.5*m2*g;

f = veiculo(L1, L2, L2eixo, m1, m2, uIz, R, g, veld, F1, F2);

xi = 0;
xf = 60;
yi = [0; 0.4; 0; -0.1];
h = 0.01;

% Método de Runge-Kutta 4 Ordem
[t_rk4, theta_rk4] = runge_kutta_4(f, h, yi, xi, xf);

figure(200 + 10 + 4);
plot(t_rk4, theta_rk4)
grid
str = {'$\theta$1','$\dot{\theta}$1','$\theta$2','$\dot{\theta}$2'};
legend(str,'Interpreter','latex');
xlabel("tempo(s)");
ylabel("posicao(rad), velocidade(rad/s)");
title(["Solucao do veiculo para m2 = 1000kg utilizando o Método de Runge-Kutta de 4 ordem e passo : 0.01"]);

