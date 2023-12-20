clear all
close all
clc

Ts = 1/1000;
y(1) = 0;
r(1) = 0;
e(1) = 0;
u(1) = 0;

% Uma solução manual
Kp = 1;
Ki = 10;

% Uma solução otimizada
x = [2.82689           20];

Kp = x(1);
Ki = x(2);



for k=2:1000

  t(k) = k*Ts;        % Tempo

  % Referência
  r(k) = 10;

  % Erro de rastreamento
  e(k) = r(k-1) - y(k-1);

  % Ação de controle
  u(k) = u(k-1)+Kp*(e(k)-e(k-1)) + Ki*Ts*e(k) ; % PID

  % discretização usando Euler forwward
  y(k) = ( y(k-1) + 357.3*Ts*u(k) ) / ( 1+28.38*Ts );

end

fonte = 21;

figure
plot(t,r,':k','LineWidth',3)
hold
plot(t,y,'b','LineWidth',3)
grid on;
xlabel('Tempo (s)',"fontsize", fonte);
ylabel('Velocidade angular (rad/s)',"fontsize", fonte);
legend('r', 'y');
set(gcf,'color','white');
h=get(gcf, "currentaxes");
set(h, "fontsize", fonte);

figure
plot(t,u,'b','LineWidth',3)
grid on;
xlabel('Tempo (s)',"fontsize", fonte);
ylabel('Tensão (V)',"fontsize", fonte);
legend('u');
set(gcf,'color','white');
h=get(gcf, "currentaxes");
set(h, "fontsize", fonte);
