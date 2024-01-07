function o = LLCL_for_opt(x)
t = zeros(1, 60000);
r = zeros(1, 60000);
d = zeros(1, 60000);
d_cos = zeros(1, 60000); 
d_sin = zeros(1, 60000); 
aux = zeros(3, 60000);
y = zeros(1, 60000); 
e0 = zeros(1, 60000); 

x_local(1) = x(1);
x_local(2) = x(2);
x_local(3) = x(3);
x_local(4) = x(4);
x_local(5) = x(5);
x_local(6) = x(6);

% Par�metros do LLCL
Lf = 1.05e-6;
Cf = 60e-6;
L1 = 200e-6;
L2 = 40e-6;

Rd = 0; % sem amortecimento passivo
% Rd = 1; % com amortecimento passivo

% sem resist�ncias parasitas
% R1 = 0;
% R2 = 0;
% com resist�ncias parasitas
R1 = 50e-3;
R2 = 50e-3;

Ts = 1/20000;

% FT em s
Gs = tf([ (Lf*Cf) (Cf*Rd) (1)], [(Cf*(Lf*(L1+L2)+(L2*L1))) (Cf*((Rd*(L1+L2))+(Lf*(R1+R2))+(L1*R2+L2*R1))) (L1+L2+(Cf*(Rd*(R1+R2)))+(R1*R2)) (R1+R2)]);

%figure
%bode(Gs)

% FT discreta
Gz = c2d(Gs,Ts,'zoh') ;

%       0.1695 z^2 - 0.08044 z - 0.03622
% Gz =  ------------------------------------
%       z^3 - 1.759 z^2 + 0.7721 z - 0.01275

[numZ, denZ] = tfdata(Gz, 'v');
b0 = numZ(2);
b1 = numZ(3);
b2 = numZ(4);
a0 = denZ(1);
a1 = denZ(2);
a2 = denZ(3);
a3 = denZ(4);

%figure
%bode(Gz)

%           b0 z^2 + b1 z + b2
% Gz =  ----------------------------
%       a0 z^3 + a1 z^2 + a2 z + a3


% Condi��es iniciais
y(1) = 0;
y(2) = 0;
y(3) = 0;
u(1) = 0;
u(2) = 0;
u(3) = 0;
e0(1) = 0;
e0(2) = 0;
e0(3) = 0;

% Dist�bio da rede
Fd = [ 0; 0; 0.3e-3 ];

% RMRAC (para cancelar o dist�rbio)
sigma_zero=0.1;
Mo=5;
m2(4)=4;
delta0=0.7;
delta1=1;
% Ganhos e par�metros para converg�ncia
% Gamma = 10;
% kappa = 2000;

% inicializa��o dos ganhos (theta) do RMRAC
% SE N�O OTIMIZAR THETAS (RODA)
ini = -1;
% Thetac(1) = ini;
% Thetas(1) = ini;
% Thetac(2) = ini;
% Thetas(2) = ini;
% Thetac(3) = ini;
% Thetas(3) = ini;
% Thetac(4) = ini;
% Thetas(4) = ini;
% Theta(:,4)= [Thetac(4); Thetas(4)];

% Controlador PI
% Kp = 1;
% Ki = 10;
upi(1) = 0;
upi(2) = 0;
upi(3) = 0;


Thetac = zeros(1, 60000);
Thetas = zeros(1, 60000);

x_local = [ x(1) x(2)  x(3)  x(4)  x(5)  x(6) ];

% Solu��es do GWO
Kp = x_local(1);
Ki = x_local(2);
Gamma = x_local(3);
kappa = x_local(4);
Thetac(4) = x_local(5);
Thetas(4) = x_local(6);

Theta(:,4)= [Thetac(4); Thetas(4)];


% inicializa��o dos ganhos (theta) do RMRAC
% ini_c = x(5);
% ini_s = x(6);
% Thetac(1) = ini_c;
% Thetas(1) = ini_s;
% Thetac(2) = ini_c;
% Thetas(2) = ini_s;
% Thetac(3) = ini_c;
% Thetas(3) = ini_s;
% Thetac(4) = ini_c;
% Thetas(4) = ini_s;
% Theta(:,4)=[ ini_c; ini_s];



for k = 4:60000
    % tempo
    t(k) = k*Ts;
    
    % corrente de refer�ncia
    if t(k)<=1
        r(k)=20*sin(2*pi*k*Ts*60);
    else
        r(k)=30*sin(2*pi*k*Ts*60);
    end
    
    %Dist�rbio ( primeira harm�nica da rede - fundamental - 60 Hz )
    d(k)=311*sin(2*pi*k*60*Ts);
    
    % rejei��o de dist�rbios
    % disturbio em fase e quadratura
    d_cos(k)=179.6*cos(2*pi*60*k*Ts);
    d_sin(k)=179.6*sin(2*pi*60*k*Ts);
    
    % sa�da do sistema (ig - corrente injetada na rede)
    if ( t(k)>=2 )
        % Dist�bio da rede
        Lg = 1e-3; % indut�ncia da rede el�trica (com varia��o param�trica)
        Fd = [ 0; 0; -1/Lg ];
    end
    
    aux(1:3,k) =  Ts*Fd*d(k);
    y(k) = (1/a0)*( -a1*y(k-1)-a2*y(k-2)-a3*y(k-3)+b0*u(k-1)+b1*u(k-2)+b2*u(k-3) ) + aux(3,k);
    
    % rejei��o usando rmrac
    e0(k) = r(k)-y(k);
    
    % Sigma modification
    norm_Theta(k)=norm(Theta(:,k));
    
    if (norm_Theta(k) < Mo)
        f_sigma_mod(k)=0;
    elseif (norm_Theta(k) >=Mo && norm_Theta(k)<= 2*Mo)
        f_sigma_mod(k)=sigma_zero*((norm_Theta(k)/Mo)-1);
    elseif (norm_Theta(k) > 2*Mo)
        f_sigma_mod(k)=sigma_zero;
    end
    
    % Lei de controle RMRAC
    ud(k)=-(Thetac(k)*d_cos(k)+Thetas(k)*d_sin(k));
    
    
    % Modulos e funcao normalizacao
    if ud(k) <= 0
        mod_u(k) = -ud(k);
    else
        mod_u(k) = ud(k);
    end
    
    if y(k) <= 0
        mod_y(k) = -y(k);
    else
        mod_y(k) = y(k);
    end
    
    % Vetor auxiliar
    Omega(:,k) = [ d_cos(k) d_sin(k) ];
    
    % Sinal majorante
    m2(k) = (1-delta0*Ts)*m2(k-1) + delta1*Ts*( 1 + mod_y(k) + mod_u(k) );
    m(k)  = m2(k)*m2(k) + (Omega(:,k)'*Gamma*Omega(:,k));
    
    % Ganhos do RMRAC
    Thetac(k+1)=Thetac(k)*(1-Ts*f_sigma_mod(k)*Gamma)-(Ts*Gamma*kappa*d_cos(k)*e0(k))/([m(k)]);
    Thetas(k+1)=Thetas(k)*(1-Ts*f_sigma_mod(k)*Gamma)-(Ts*Gamma*kappa*d_sin(k)*e0(k))/([m(k)]);
    
    Theta(:,k+1) = [ Thetac(k+1); Thetas(k+1) ];
    
    % A��o de controle do controlador PI
    upi(k) = upi(k-1)+Kp*( e0(k)-e0(k-1) ) + Ki*Ts*e0(k);
    
    % A��o de controle total
    u(k) = upi(k) + ud(k);
    
    if isnan(e0)
        e0 = 1000;
    end
    
end

o = mean(abs(e0));

end