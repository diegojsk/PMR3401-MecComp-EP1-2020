%% PMR 3401 - Mec�nica Computacional para Mecatr�nica - 2019
%  Exerc�cio Programa 1
%
%  Autores: Diego Jun Sato Kurashima
%           Felipe Gomes de Melo D'Elia
%
%% Teste
%
% y' = x => y = (x^2)/2

f = @(x,y) x;
euler_method(f,1,0,0,5)
runge_kutta_2(f,1,0,0,5)
runge_kutta_4(f,1,0,0,5)

%% Equa��o do Modelo
%
% Setar par�metros
L1 = 2;
L2 = 2.5;
L2eixo = 1.8;
m1 = 450;
m2 = 650;
F1 = -0.5;
F2 = -0.5;
uIz = 2.7;
R = 0.3;
g = 9.81;
veld = 80;

% Equa��es 
%
% thetha1 = y(1)
% thetha1' = y(2)
%
% thetha2 = y(3)
% thetha2' = y(4)

% Equa��o de Movimento
f = @(t,y) [y(2);
            y(2)*2;
            y(3);
            y(4) + y(1)];
        
euler_method(f,1,[0; 0; 0; 0],0,5)
