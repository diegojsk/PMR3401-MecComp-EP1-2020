%% Teste
% Testar Execução dos Métodos Numéricos: 
% y' = 2x => y = x^2

t = @(x,y) 2*x;
euler_method(t,1,1,1,5)
runge_kutta_2(t,1,1,1,5)
runge_kutta_4(t,1,1,1,5)