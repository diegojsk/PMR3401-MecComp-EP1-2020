%% PMR 3401 - Mecânica Computacional para Mecatrônica - 2020
%  Exercício Programa 1
%  Script para analisar solucoes referentes a Questao 1
%
%  Autores: Diego Jun Sato Kurashima - 10274231
%           Felipe Gomes de Melo D'Elia - 10340624
%
clear all
clc

%% Questao 1

% Setar parâmetros
L1 = 2;
L2 = 2.5;
L2eixo = 1.8;
m1 = 450;
m2 = 650;
uIz = 2.7;
R = 0.3;
g = 9.81;
veld = 80/3.6;
F1 = -0.5*m1*g;
F2 = -0.5*m2*g;

f = @(t,y) [ y(2);
            (1/((L1^2)*L2*R*(m2*cos(2*y(1) - 2*y(3)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*y(1) - 2*y(3)))*(y(2)^2) + (2*L1*(L2^2)*R*m2*sin(y(1)-y(3)))*(y(4)^2) + (-2*L2*uIz*veld)*y(2) + (-2*L1*uIz*veld*cos(y(1) - y(3)))*y(4) + (-R*L1*( L2eixo*F2*sin(y(1) - 2*y(3)) + 2*sin(y(1))*( F1*L2 + (1/2)*L2eixo*F2))));
             y(4);
          ( (1/((L2^2)*R*m2))*( (-L1*L2*R*m2*cos(y(1) - y(3)) )*( (1/((L1^2)*L2*R*(m2*cos(2*y(1) - 2*y(3)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*y(1) - 2*y(3)))*(y(2)^2) + (2*L1*(L2^2)*R*m2*sin(y(1)-y(3)))*(y(4)^2) + (-2*L2*uIz*veld)*y(2) + (-2*L1*uIz*veld*cos(y(1) - y(3)))*y(4) + (-R*L1*( L2eixo*F2*sin(y(1) - 2*y(3)) + 2*sin(y(1))*( F1*L2 + (1/2)*L2eixo*F2))))) + (L1*L2*R*m2*sin(y(1) - y(3)))*(y(2)^2) + (-uIz*veld)*y(4) + (L2eixo*sin(y(3))*R*F2) ))];

% Condicoes Iniciais

xi = 0;
xf = 60;
yi = [0; 0.4; 0; -0.1];

% Testar a resolução para diferentes valores de passo h

i = 1; % variável para identificacao dos graficos

% Analisando para vários passos h
for h = [0.5, 0.1, 0.05, 0.01]

    % Método de Euler
    [t_e, theta_e] = euler_method(f, h, yi, xi, xf);

    figure(100 + i*10 + 1);
    plot(t_e, theta_e)
    grid
    hold on

    % Calcular as aceleracoes
    for n = 1:((xf-xi)/h)+1
        theta_e_ll(:,n) = [(1/((L1^2)*L2*R*(m2*cos(2*theta_e(1,n) - 2*theta_e(3,n)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*theta_e(1,n) - 2*theta_e(3,n)))*(theta_e(2,n)^2) + (2*L1*(L2^2)*R*m2*sin(theta_e(1,n)-theta_e(3,n)))*(theta_e(4,n)^2) + (-2*L2*uIz*veld)*theta_e(2,n) + (-2*L1*uIz*veld*cos(theta_e(1) - theta_e(3,n)))*theta_e(4,n) + (-R*L1*( L2eixo*F2*sin(theta_e(1,n) - 2*theta_e(3,n)) + 2*sin(theta_e(1,n))*( F1*L2 + (1/2)*L2eixo*F2))));
                          ( (1/((L2^2)*R*m2))*( (-L1*L2*R*m2*cos(theta_e(1,n) - theta_e(3,n)) )*( (1/((L1^2)*L2*R*(m2*cos(2*theta_e(1,n) - 2*theta_e(3,n)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*theta_e(1,n) - 2*theta_e(3,n)))*(theta_e(2,n)^2) + (2*L1*(L2^2)*R*m2*sin(theta_e(1,n)-theta_e(3,n)))*(theta_e(4,n)^2) + (-2*L2*uIz*veld)*theta_e(2,n) + (-2*L1*uIz*veld*cos(theta_e(1) - theta_e(3,n)))*theta_e(4,n) + (-R*L1*( L2eixo*F2*sin(theta_e(1,n) - 2*theta_e(3,n)) + 2*sin(theta_e(1,n))*( F1*L2 + (1/2)*L2eixo*F2))))) + (L1*L2*R*m2*sin(theta_e(1,n) - theta_e(3,n)))*(theta_e(2,n)^2) + (-uIz*veld)*theta_e(4,n) + (L2eixo*sin(theta_e(3,n))*R*F2) ))];
    end
    % Gráfico completo
    plot(t_e, theta_e_ll);
    grid
    str = {'$\theta$1','$\dot{\theta}$1','$\theta$2','$\dot{\theta}$2','$\ddot{\theta}$1','$\ddot{\theta}$2'};
    legend(str,'Interpreter','latex');
    xlabel("tempo(s)");
    ylabel("posicao(rad), velocidade(rad/s) ou aceleracao(rad/s^2)");
    title(["Solucao completa do veiculo utilizando o Método de Euler e passo :" h]);
    hold off


    % Método de Runge-Kutta 2 Ordem
    [t_rk2, theta_rk2] = runge_kutta_2(f, h, yi, xi, xf);

    figure(100 + i*10 + 2);
    plot(t_rk2, theta_rk2)
    grid
    hold on

    for n = 1:((xf-xi)/h)+1
        % Calcular as aceleracoes
        theta_rk2_ll(:,n) = [(1/((L1^2)*L2*R*(m2*cos(2*theta_rk2(1,n) - 2*theta_rk2(3,n)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*theta_rk2(1,n) - 2*theta_rk2(3,n)))*(theta_rk2(2,n)^2) + (2*L1*(L2^2)*R*m2*sin(theta_rk2(1,n)-theta_rk2(3,n)))*(theta_rk2(4,n)^2) + (-2*L2*uIz*veld)*theta_rk2(2,n) + (-2*L1*uIz*veld*cos(theta_rk2(1) - theta_rk2(3,n)))*theta_rk2(4,n) + (-R*L1*( L2eixo*F2*sin(theta_rk2(1,n) - 2*theta_rk2(3,n)) + 2*sin(theta_rk2(1,n))*( F1*L2 + (1/2)*L2eixo*F2))));
                          ( (1/((L2^2)*R*m2))*( (-L1*L2*R*m2*cos(theta_rk2(1,n) - theta_rk2(3,n)) )*( (1/((L1^2)*L2*R*(m2*cos(2*theta_rk2(1,n) - 2*theta_rk2(3,n)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*theta_rk2(1,n) - 2*theta_rk2(3,n)))*(theta_rk2(2,n)^2) + (2*L1*(L2^2)*R*m2*sin(theta_rk2(1,n)-theta_rk2(3,n)))*(theta_rk2(4,n)^2) + (-2*L2*uIz*veld)*theta_rk2(2,n) + (-2*L1*uIz*veld*cos(theta_rk2(1) - theta_rk2(3,n)))*theta_rk2(4,n) + (-R*L1*( L2eixo*F2*sin(theta_rk2(1,n) - 2*theta_rk2(3,n)) + 2*sin(theta_rk2(1,n))*( F1*L2 + (1/2)*L2eixo*F2))))) + (L1*L2*R*m2*sin(theta_rk2(1,n) - theta_rk2(3,n)))*(theta_rk2(2,n)^2) + (-uIz*veld)*theta_rk2(4,n) + (L2eixo*sin(theta_rk2(3,n))*R*F2) ))];
    end
    plot(t_rk2, theta_rk2_ll)
    grid
    str = {'$\theta$1','$\dot{\theta}$1','$\theta$2','$\dot{\theta}$2','$\ddot{\theta}$1','$\ddot{\theta}$2'};
    legend(str,'Interpreter','latex');
    xlabel("tempo(s)");
    ylabel("posicao(rad), velocidade(rad/s) ou aceleracao(rad/s^2)");
    title(["Solucao completa do veiculo utilizando o Método de Runge-Kutta de 2 ordem e passo :" h]);
    hold off


    % Método de Runge-Kutta 4 Ordem
    [t_rk4, theta_rk4] = runge_kutta_4(f, h, yi, xi, xf);

    figure(100 + i*10 + 4);
    plot(t_rk4, theta_rk4)
    grid
    hold on
    for n = 1:((xf-xi)/h)+1
        % Calcular as aceleracoes
        theta_rk4_ll(:,n) = [(1/((L1^2)*L2*R*(m2*cos(2*theta_rk4(1,n) - 2*theta_rk4(3,n)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*theta_rk4(1,n) - 2*theta_rk4(3,n)))*(theta_rk4(2,n)^2) + (2*L1*(L2^2)*R*m2*sin(theta_rk4(1,n)-theta_rk4(3,n)))*(theta_rk4(4,n)^2) + (-2*L2*uIz*veld)*theta_rk4(2,n) + (-2*L1*uIz*veld*cos(theta_rk4(1) - theta_rk4(3,n)))*theta_rk4(4,n) + (-R*L1*( L2eixo*F2*sin(theta_rk4(1,n) - 2*theta_rk4(3,n)) + 2*sin(theta_rk4(1,n))*( F1*L2 + (1/2)*L2eixo*F2))));
                          ( (1/((L2^2)*R*m2))*( (-L1*L2*R*m2*cos(theta_rk4(1,n) - theta_rk4(3,n)) )*( (1/((L1^2)*L2*R*(m2*cos(2*theta_rk4(1,n) - 2*theta_rk4(3,n)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*theta_rk4(1,n) - 2*theta_rk4(3,n)))*(theta_rk4(2,n)^2) + (2*L1*(L2^2)*R*m2*sin(theta_rk4(1,n)-theta_rk4(3,n)))*(theta_rk4(4,n)^2) + (-2*L2*uIz*veld)*theta_rk4(2,n) + (-2*L1*uIz*veld*cos(theta_rk4(1) - theta_rk4(3,n)))*theta_rk4(4,n) + (-R*L1*( L2eixo*F2*sin(theta_rk4(1,n) - 2*theta_rk4(3,n)) + 2*sin(theta_rk4(1,n))*( F1*L2 + (1/2)*L2eixo*F2))))) + (L1*L2*R*m2*sin(theta_rk4(1,n) - theta_rk4(3,n)))*(theta_rk4(2,n)^2) + (-uIz*veld)*theta_rk4(4,n) + (L2eixo*sin(theta_rk4(3,n))*R*F2) ))];
    end

    plot(t_rk4, theta_rk4_ll)
    grid
    str = {'$\theta$1','$\dot{\theta}$1','$\theta$2','$\dot{\theta}$2','$\ddot{\theta}$1','$\ddot{\theta}$2'};
    legend(str,'Interpreter','latex');
    xlabel("tempo(s)");
    ylabel("posicao(rad), velocidade(rad/s) ou aceleracao(rad/s^2)");
    title(["Solucao completa do veiculo utilizando o Método de Runge-Kutta de 4 ordem e passo :" h]);
    hold off

    i = i + 1;
end
