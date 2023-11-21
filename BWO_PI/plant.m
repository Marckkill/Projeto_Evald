function o = plant(x)

Kp = x(1);
Ki = x(2);


% Taxa de amostragem
Ts = 1/1000;
% Condições iniciais (pq o for começa em k=2
y(1) = 0;
r(1) = 0;
e(1) = 0;
u(1) = 0;

for k=2:1002

  % Tempo
  t(k) = k*Ts;        

  % Referência
  r(k) = 10;

  % Erro de rastreamento
  e(k) = r(k-1) - y(k-1);

  % Ação de controle
  u(k) = u(k-1)+Kp*(e(k)-e(k-1)) + Ki*Ts*e(k) ; % PID

  % discretização usando Euler forwward
  y(k) = ( y(k-1) + 357.3*Ts*u(k) ) / ( 1+28.38*Ts );

end

o = mean(abs(e));

end