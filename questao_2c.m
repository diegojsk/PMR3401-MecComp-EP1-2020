%% PMR 3401 - Mecânica Computacional para Mecatrônica - 2020
%  Exercício Programa 1
%  Script para analisar solucoes referentes a Questao 2c 
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 
%
clear all
clc
%% Questao 2a

L1 = 2;
L2 = 2.5;
L2eixo = 1.8;
m1 = 450;
m2 = 650;
uIz = 2.7;
R = 0.3;
g = 9.81;
veld = 120/3.6;
F1 = -0.5*m1*g;
F2 = -0.5*m2*g;

f = @(t,y) [ y(2);
            (1/((L1^2)*L2*R*(m2*cos(2*y(1) - 2*y(3)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*y(1) - 2*y(3)))*(y(2)^2) + (2*L1*(L2^2)*R*m2*sin(y(1)-y(3)))*(y(4)^2) + (-2*L2*uIz*veld)*y(2) + (-2*L1*uIz*veld*cos(y(1) - y(3)))*y(4) + (-R*L1*( L2eixo*F2*sin(y(1) - 2*y(3)) + 2*sin(y(1))*( F1*L2 + (1/2)*L2eixo*F2))));
             y(4);
          ( (1/((L2^2)*R*m2))*( (-L1*L2*R*m2*cos(y(1) - y(3)) )*( (1/((L1^2)*L2*R*(m2*cos(2*y(1) - 2*y(3)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*y(1) - 2*y(3)))*(y(2)^2) + (2*L1*(L2^2)*R*m2*sin(y(1)-y(3)))*(y(4)^2) + (-2*L2*uIz*veld)*y(2) + (-2*L1*uIz*veld*cos(y(1) - y(3)))*y(4) + (-R*L1*( L2eixo*F2*sin(y(1) - 2*y(3)) + 2*sin(y(1))*( F1*L2 + (1/2)*L2eixo*F2))))) + (L1*L2*R*m2*sin(y(1) - y(3)))*(y(2)^2) + (-uIz*veld)*y(4) + (L2eixo*sin(y(3))*R*F2) ))];

      
xi = 0;
xf = 60;
yi = [0; 0.4; 0; -0.1];
h = 0.01;

% Método de Runge-Kutta 4 Ordem
[t_rk4, theta_rk4] = runge_kutta_4(f, h, yi, xi, xf);

figure(200 + 30 + 4);
plot(t_rk4, theta_rk4)
grid
str = {'$\theta$1','$\dot{\theta}$1','$\theta$2','$\dot{\theta}$2'};
legend(str,'Interpreter','latex');
xlabel("tempo(s)");
ylabel("posicao(rad), velocidade(rad/s)");
title(["Solucao do veiculo para Vel_d = 120km/h utilizando o Método de Runge-Kutta de 4 ordem e passo : 0.01"]);     
