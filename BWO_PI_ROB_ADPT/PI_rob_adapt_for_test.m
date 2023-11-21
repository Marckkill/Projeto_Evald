close all
clear all
clc

format long g

%% Parte 1 - Parametros de Simulacao e Definicao Plantas
fs=5.04e3;       % Frequencia de discretizacao ( = Passo da simulacao)
Ts=1/fs;         % Periodo de Amostragem (s)
Ttotal=3;        % Tempo total de simulacao (s)
tempo=0:Ts:Ttotal;

Vlink = 500;   % tensao de barramento (monofasico 127V rede)
Cf = 62e-6;    % capacitor filtro LCL

rc = 50e-3;    % resistencia referente ao primeiro braco do LCL
rg = 50e-3;    % resistencia referente segundo braco do LCL

%Rd = 1;        % resistencia de amortecimento passivo (serie capacitor LCL)
Rd = 0;
Ro = 0;        % curto-circuito terminais fonte (sem considerar fonte Vg)

Lc = 1e-3;     % Indutancia referente ao primeiro braco do LCL
Lg = 0.3e-3;   % Indutancia referente ao segundo braco do LCL

% Modelo de terceira ordem - para controle da corrente do lado da rede (ig)
% - conforme a Tese - Anexo C - Modelagem da planta
s = tf('s');
sysC3 = 1/(rc + rg + Lc*s + Lg*s + Cf*Lc*Lg*s^3 + Cf*Lc*rg*s^2 + ...
        Cf*Lg*rc*s^2 + Cf*rc*rg*s);
[numC3,denC3] = tfdata(sysC3,'v');

Lg_grid = 1e-3 ; %10e-3;       % indutancia da rede
% Modelo da planta completa considerando a variacao parametrica
sysC3_VP = 1/(rc + rg + Lc*s + (Lg+Lg_grid)*s + Cf*Lc*(Lg+Lg_grid)*s^3 + Cf*Lc*rg*s^2 + ...
        Cf*(Lg+Lg_grid)*rc*s^2 + Cf*rc*rg*s);
[numC3_VP,denC3_VP] = tfdata(sysC3_VP,'v');    

% Modelo considerando Rd
sysC3_Rd = ((Cf * Rd * s + 1) ) / (rc + rg + Ro + Lc * s ...
        + Lg * s + Cf * Lc * Lg * s^3 + Cf * Lc * rg * s^2 + Cf * Lg ...
         * rc * s^2 + Cf * Lc * Rd * s^2 + Cf * Lg * Rd * s^2 + ...
         Cf * Lc * Ro * s^2 + Cf * rc * rg * s + Cf * rc * rg * s ...
         + Cf * rg * Rd * s + Cf * rc * Ro * s + Cf * Rd * Ro * s);

% Gp(s) - Parte Modelada da Planta - 1 ordem
% negligenciando os polos complexos conjugados
kpc = 1 / (Lc + Lg);
numC = kpc;
denC = [1 (rg + rc) / (Lc + Lg) ];
sysC = tf( numC , denC);

% Polo real do modelo de primeira ordem
p4 = abs( roots (denC) );

% Transformada Z, com ZOH, aplicada a planta Gp(s) (reduzido) 1 ordem
Gs1 = tf (numC , denC);
Gz1 = c2d (Gs1 , Ts);
[numZ , denZ] = tfdata( Gz1 , 'v');
b0 = numZ(1);
b1 = numZ(2);
a0 = denZ(1);
a1 = denZ(2);

% Transformada Z, com ZOH, aplicada a planta Gp(s) (completo) 3.a ordem
Gs2=tf(numC3,denC3);
Gz2=c2d(Gs2,Ts);
[numZ3,denZ3]=tfdata(Gz2,'v');
d0=numZ3(1);
d1=numZ3(2);
d2=numZ3(3);
d3=numZ3(4);
c0=denZ3(1);
c1=denZ3(2);
c2=denZ3(3);
c3=denZ3(4);

% Transformada Z, com ZOH, aplicada a planta Gp(s) (completo) 3.a ordem
% com variacao parametrica
Gs3=tf(numC3_VP,denC3_VP);
Gz3=c2d(Gs3,Ts);
[numZ3,denZ3]=tfdata(Gz3,'v');
d0_VP=numZ3(1);
d1_VP=numZ3(2);
d2_VP=numZ3(3);
d3_VP=numZ3(4);
c0_VP=denZ3(1);
c1_VP=denZ3(2);
c2_VP=denZ3(3);
c3_VP=denZ3(4);

delta0=0.7;
delta1=1;


%% Parte 2 - Parametros PI-RMRAC (sem Wm)

% Parametros da lei de adaptacao parametrica
% Inicializacao dos vetores de entrada, saida e sinais internos (memoria)
y=zeros(1,length(tempo));           % Saida desejada
y1=zeros(1,length(tempo));          % Saida do modelo da planta - 1 ordem
y3=zeros(1,length(tempo));          % Saida do modelo da planta - 3 ordem
u=zeros(1,length(tempo));           % Lei de controle
r=zeros(1,length(tempo));           % Sinal de referencia
e0=zeros(1,length(tempo));          % Erro de regulação
m=zeros(1,length(tempo));           % Normalizador m
m2=4*ones(1,length(tempo));         % Normalizador m
mod_y=ones(1,length(tempo));        % Normalizador m
mod_u=ones(1,length(tempo));        % Normalizador m
vg=zeros(1,length(tempo));          % Tensao da rede (disturbio periodico)
d_cos=zeros(1,length(tempo));       % Decomposicao em fase e quadradatura
d_sin=zeros(1,length(tempo));       % do disturbio periodico (sin e cos)
norm_Theta=zeros(1,length(tempo));
sigma=zeros(1,length(tempo));
sigma_zero=0.1;
Mo=5;

% Valor Inicial dos Ganhos de adaptacao
MAO = [ 250.5551      833.5247      17.16772     -8.244002     0.8741073     -10.09149     -1.854357     -8.684884 ];
JSOA = [ 1.47361      500.8178     0.9126044    -0.3160528    -0.6875876      -1.23594      -0.05313    -0.5745087 ];

x = MAO;
x = JSOA;

Gamma = x(1);
kappa = x(2);
Theta1 = x(3)*ones(1,length(tempo));
Theta2 = x(4)*ones(1,length(tempo));
Theta3 = x(5)*ones(1,length(tempo));
Theta4 = x(6)*ones(1,length(tempo));
Thetac = x(7)*ones(1,length(tempo));
Thetas = x(8)*ones(1,length(tempo));


Theta=[Theta1; Theta2; Theta3; Theta4; Thetac; Thetas];


%% Parte 3 - Simulacao do Sistema em Malha Fechada (PI-RMRAC)

for k=4:length(tempo)
        
    % Sinal de Referencia 
    if tempo(k) < 1
        r(k)=10*sin(2*pi*60*Ts*k);
    elseif tempo(k) < 1.75
        r(k)=20*sin(2*pi*60*Ts*k);
    else
        r(k)=30*sin(2*pi*60*Ts*k);
    end
    
    % Definicao do sinal de Tensao da Rede e sinais de decomposicao em fase
    % e quadratura -> d_sin(k) e d_cos(k)
    % phi=pi/4;
    phi = 0;
    vg(k)=179.6*sin(2*pi*60*k*Ts+phi);
    d_cos(k)=179.6*cos(2*pi*60*k*Ts);
    d_sin(k)=179.6*sin(2*pi*60*k*Ts);
    
    % Saida da planta do modelo completo de 1.a ordem - Somente Teste
    y1(k)=-y1(k-1)*(a1)+u(k-1)*(b1)-b1*vg(k-1);
    
    % Saida da planta do modelo completo de 3.a ordem - Planta Real
    if tempo(k) < 2.5
        y3(k)=(1/c0)*(-c3*y3(k-3)-c2*y3(k-2)-c1*y3(k-1)+ ...
            d3*u(k-3)+d2*u(k-2)+d1*u(k-1)+d0*u(k))-(d3*vg(k-1)+d2*vg(k-1)+ ...
            d1*vg(k-1)+d0*vg(k-1));        
    else
        y3(k)=(1/c0_VP)*(-c3_VP*y3(k-3)-c2_VP*y3(k-2)-c1_VP*y3(k-1)+ ...
            d3_VP*u(k-3)+d2_VP*u(k-2)+d1_VP*u(k-1)+d0_VP*u(k))-(d3_VP*vg(k-1)+d2_VP*vg(k-1)+ ...
            d1_VP*vg(k-1)+d0_VP*vg(k-1));        
    end
    
    % Escolha do modelo da planta (descomentar o utilizado)
%   y(k)=y1(k);
    y(k)=y3(k);
     
    % Erro de regulação (e0)
    e0(k)=r(k)-y(k);
    
    % Sigma modification
    norm_Theta(k)=norm(Theta(:,k));

    if (norm_Theta(k) < Mo)
        sigma(k)=0;
    elseif (norm_Theta(k) >=Mo && norm_Theta(k)<= 2*Mo)
        sigma(k)=sigma_zero*((norm_Theta(k)/Mo)-1);
    elseif (norm_Theta > 2*Mo)
        sigma(k)=sigma_zero;
    end
       

    % Lei de controle utilizada (descomentar a que for de interesse)
    u(k)=-(Theta2(k)*u(k-1)+Theta3(k)*y(k)+Theta4(k)*e0(k-1)+r(k))/Theta1(k)-(Thetac(k)*d_cos(k)+Thetas(k)*d_sin(k))/Theta1(k);

    % Modulos e funcao normalizacao
    if u(k) < 0
        mod_u(k) = -u(k);
    else
        mod_u(k) = u(k);        
    end

    if y(k) < 0
        mod_y(k) = -y(k);
    else
        mod_y(k) = y(k);    
    end

    Omega(:,k) = [  u(k) u(k-1) y(k) e0(k-1) d_cos(k) d_sin(k) ];
    
    m2(k) = (1-delta0*Ts)*m2(k-1) + delta1*Ts*( 1 + mod_y(k) + mod_u(k) );
    m(k)  = m2(k)*m2(k) + (Omega(:,k)'*Gamma*Omega(:,k));
    
    % Ganhos adaptativos
    Theta1(k+1)=Theta1(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*u(k)*e0(k))/([m(k)]);
    Theta2(k+1)=Theta2(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*u(k-1)*e0(k))/([m(k)]);
    Theta3(k+1)=Theta3(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*y(k)*e0(k))/([m(k)]);
    Theta4(k+1)=Theta4(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*e0(k-1)*e0(k))/([m(k)]);
    Thetac(k+1)=Thetac(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*d_cos(k)*e0(k))/([m(k)]);
    Thetas(k+1)=Thetas(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*d_sin(k)*e0(k))/([m(k)]);
    Theta(:,k+1)=[Theta1(k+1);Theta2(k+1);Theta3(k+1);Theta4(k+1);Thetac(k+1);Thetas(k+1)]; 
    
end
 
%% Parte 4 - Graficos 

fonte = 21;

a=figure;
% plot(tempo,r,'k','LineWidth',2);
plot(tempo,r,':k','LineWidth',4);
hold on
% plot(t,y,'color',[0,0,0]+0.7,'LineWidth',2);
plot(tempo,y,'b','LineWidth',2); %4
grid on;
h2 = legend('$r$','$y$',([345, 100, 0, 0]));
% h2 = legend('$y_{m}$','$y$',([345, 100, 0, 0]));
set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
set(gcf,'Units','centimeters','PaperSize',[13 7]);   
% Recortar a figura da página
set(gca,'fontsize',fonte,'units','norm');
ylabel('Current (A)','fontsize',fonte);
xlabel('Time (s)','fontsize',fonte);
% xlim([0 0.030])
% ylim([-25 25])
% xlim([0.99 1.03])
% ylim([-25 25])
% xlim([1.74 1.78])
% ylim([-30 30])
% xlim([2.49 2.52])
% ylim([-35 35])


figure
plot(tempo,e0,'b','LineWidth',3);
hold
% plot(tempo,E1,'color',[0,0,0]+0.7,'LineWidth',2);
% plot(tempo,E1,':k','LineWidth',4);
grid on;
% h2 = legend('$e_{1}$','$\epsilon$',([345, 100, 0, 0]));
h2 = legend('$e_{0}$',([345, 100, 0, 0]));
set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
set(gcf,'Units','centimeters','PaperSize',[13 7]);                                     % Recortar a figura da página
set(gca,'fontsize',fonte,'units','norm');
ylabel('Errors (A)','fontsize',fonte);
xlabel('Time (s)','fontsize',fonte);
% xlim([0 0.03])
% ylim([-35 35])
% xlim([0.99 1.05])
% ylim([-2 2])
% xlim([1.74 1.82])
% ylim([-2 2])
% xlim([2.48 2.6])
% ylim([-6 6])

figure
plot(tempo,u,'b','LineWidth',2);
grid on;
h2 = legend('$u$',([345, 100, 0, 0]));
set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
set(gcf,'Units','centimeters','PaperSize',[13 7]);                                     % Recortar a figura da página
set(gca,'fontsize',fonte,'units','norm');
ylabel('Control action (V)','fontsize',fonte);
xlabel('Time (s)','fontsize',fonte);
% ylim([-0.08 0.08])
% xlim([95 130])
% ylim([-0.08 0.08])
% xlim([145 170])

tempo(k+1) = tempo(k);

figure
set(gcf,'renderer','painters')
plot(tempo,Theta(1,:),'k','LineWidth',4);
hold on
plot(tempo,Theta(2,:),':b','LineWidth',3);
plot(tempo,Theta(3,:),'--r','LineWidth',2);
plot(tempo,Theta(4,:),'-.c','LineWidth',2);
plot(tempo,Theta(5,:),'dg','LineWidth',1);
plot(tempo,Theta(6,:),'color',[0,0,0]+0.6,'LineWidth',3);
grid on;
h2 = legend('$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$','$\theta_{c}$','$\theta_{s}$',([345, 100, 0, 0]));
set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
set(gcf,'Units','centimeters','PaperSize',[14 7]);                                     % Recortar a figura da página
set(gca,'fontsize',fonte,'units','norm');
ylabel('Gains','fontsize',fonte);
xlabel('Time (s)','fontsize',fonte);
% ylim([-4 8])