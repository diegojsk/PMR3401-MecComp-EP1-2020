%% Teste
% Arquivo com testes para validar a implementa��o

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

%% Testar plotagem de gr�ficos
