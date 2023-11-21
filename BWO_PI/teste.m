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
Kd = 0.001;

% Uma solução otimizada
x = [ 2.59091           20  0.000267286 ];
x = [ 2.75783           20  0.000100363 ];
x = [ 2.74551           20  0.000112692 ] ;
x = [ 2.59786           20  0.000260333 ];
x = [ 2.73023           20  0.000127972 ];
x = [ 2.59827           20  0.000259929 ];
x = [ 2.59091           20  0.000267286 ];

Kp = x(1);
Ki = x(2);
Kd = x(3);



for k=2:1000

  t(k) = k*Ts;        % Tempo

  % Referência
  r(k) = 10;

  % Erro de rastreamento
  e(k) = r(k-1) - y(k-1);

  % Ação de controle
  u(k) = u(k-1)+Kp*(e(k)-e(k-1)) + Ki*Ts*e(k) + Kd*(e(k)-e(k-1))/Ts ; % PID

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
