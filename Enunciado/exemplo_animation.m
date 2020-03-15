% Este script gera uma anima??o do movimento dos ve?culos
% Sendo: 
%(xi,yi): A posi??o da massa i=(1,2)
% Li: Comprimento do carro i=(1,2)
% L2eixo: posi??o do eixo do segundo carro
% ? preciso que esteja declarado:
% Um vetor de tempo t e os g.d.l theta(1 e 2,t), ouseja, Uma das linhas tem theta1 e
% ourta o theta2 com tempo nas colunas.

% Este script gera uma anima??o cujo exemplo consta no moodle. Essa
% anima??o n?o precisa constar no trabalho impresso e apenas serve para
% ajudar a interpreta??o dos resultados caso haja interesse.



for ti=1:30:numel(t)
    
    figure(2);
    clf(2);
    
    x1=L1*sin(theta(1,ti));
    y1=-L1*cos(theta(1,ti));
    x2=L1*sin(theta(1,ti)) +L2*sin(theta(2,ti));
    y2=-L1*cos(theta(1,ti)) -L2*cos(theta(2,ti));
    x1e=L1*sin(theta(1,ti)) - l*cos(theta(1,ti));
    y1e=- L1*cos(theta(1,ti)) - l*sin(theta(1,ti));
    x1d=L1*sin(theta(1,ti)) + l*cos(theta(1,ti));
    y1d=- L1*cos(theta(1,ti)) + l*sin(theta(1,ti));
    
    x2e=L1*sin(theta(1,ti)) +L2eixo*sin(theta(2,ti))-cos(theta(2,ti))*l;
    y2e=-L1*cos(theta(1,ti)) -L2eixo*cos(theta(2,ti))-sin(theta(2,ti))*l;
    x2d=L1*sin(theta(1,ti)) +L2eixo*sin(theta(2,ti))+cos(theta(2,ti))*l;
    y2d=-L1*cos(theta(1,ti)) -L2eixo*cos(theta(2,ti))+sin(theta(2,ti))*l;
    
    
    plot([0 x1],[0 y1],'k-o');
    hold on;
    plot([-l l],[0 0],'k-o');
    plot([x1e x1d],[y1e y1d],'k-o');
    plot([x2e x2d],[y2e y2d],'k-o');
    plot([x1 x2],[y1 y2],'k-o');
    title(['Tempo:' num2str(t(ti)) 's']);
    
    axis([-(+1.25*(L1+L2)+0.2)/2 (+1.25*(L1+L2)+0.2)/2 -1.25*(L1+L2) 0.2]);
%     grid on
    drawnow
    pause(0.1);
%     close(2);
    
end