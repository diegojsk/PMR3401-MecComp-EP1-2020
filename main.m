%% PMR 3401 - Mec�nica Computacional para Mecatr�nica - 2020
%  Exerc�cio Programa 1
%
%  Autores: Diego Jun Sato Kurashima
%           Felipe Gomes de Melo D'Elia
%
%% Teste
% Testar Execu��o dos M�todos Num�ricos: 
% y' = 2x => y = x^2

t = @(x,y) 2*x;
euler_method(t,1,1,1,5);
runge_kutta_2(t,1,1,1,5);
runge_kutta_4(t,1,1,1,5);

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
veld = 80/3.6;

% Equa��es 
%
% thetha1 = y(1)
% thetha1' = y(2)
%
% thetha2 = y(3)
% thetha2' = y(4)

% Equa��o de Movimento

% A0 = ((L1^2)*L2*R*(m2*cos(2*y(1) - 2*y(3)) - 2*m1 - m2));
% A1 = ((L1^2)*L2*R*m2*sin(2*y(1) - 2*y(3)));
% A2 = (2*L1*(L2^2)*R*m2*sin(y(1)-y(3)));
% A3 = (-2*L2*uIz*veld);
% A4 = (-2*L1*uIz*veld*cos(y(1) - y(3)));
% A5 = (-R*L1*( L2eixo*F2*sin(y(1) - 2*y(3)) + 2*sin(y(1))*( F1*L2 + (1/2)*L2eixo*F2)));

% B0 = ((L2^2)*R*m2);
% B1 = (-L1*L2*R*m2*cos(y(1) - y(3)));
% B2 = (L1*L2*R*m2*sin(y(1) - y(3)));
% B3 = (-uIz*veld);
% B4 = (L2eixo*sin(y(3))*R*F2);


% f = @(t,y) [ y(2);
%             (1/A0)*( A1*(y(2)^2) + A2*(y(4)^2) + A3*y(2) + A4*y(2) + A5);
%              y(3);
%            ((1/B0)*( B1*((1/A0)*( A1*(y(2)^2) + A2*(y(4)^2) + A3*y(2) + A4*y(2)+ A5)) + B2*(y(2)^2) + B3*y(4) + B4 ))];
           
f = @(t,y) [ y(2);
            (1/((L1^2)*L2*R*(m2*cos(2*y(1) - 2*y(3)) - 2*m1 - m2)))*(  ((L1^2)*L2*R*m2*sin(2*y(1) - 2*y(3)))*(y(2)^2) + (2*L1*(L2^2)*R*m2*sin(y(1)-y(3)))*(y(4)^2) + (-2*L2*uIz*veld)*y(2) + (-2*L1*uIz*veld*cos(y(1) - y(3)))*y(2) + (-R*L1*( L2eixo*F2*sin(y(1) - 2*y(3)) + 2*sin(y(1))*( F1*L2 + (1/2)*L2eixo*F2))));
             y(3);
           ((1/((L2^2)*R*m2))*( (-L1*L2*R*m2*cos(y(1) - y(3)))*((1/((L1^2)*L2*R*(m2*cos(2*y(1) - 2*y(3)) - 2*m1 - m2)))*( ((L1^2)*L2*R*m2*sin(2*y(1) - 2*y(3)))*(y(2)^2) + (2*L1*(L2^2)*R*m2*sin(y(1)-y(3)))*(y(4)^2) + (-2*L2*uIz*veld)*y(2) + (-2*L1*uIz*veld*cos(y(1) - y(3)))*y(2) + (-R*L1*( L2eixo*F2*sin(y(1) - 2*y(3)) + 2*sin(y(1))*( F1*L2 + (1/2)*L2eixo*F2))))) + (L1*L2*R*m2*sin(y(1) - y(3)))*(y(2)^2) + (-uIz*veld)*y(4) + (L2eixo*sin(y(3))*R*F2) ))];
           
        
euler_method(f,1,[0; 0.4; 0; -0.5],0,60)
runge_kutta_2(f,1,[0; 0.4; 0; -0.5],0,60)
runge_kutta_4(f,1,[0; 0.4; 0; -0.5],0,60)
