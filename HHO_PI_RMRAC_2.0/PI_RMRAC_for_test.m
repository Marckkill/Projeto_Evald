close all
clear all
clc

format long g

%% Parte 1 - Parametros de Simulacao e Definicao Plantas
fs=5.04e3;       % Frequencia de discretizacao ( = Passo da simulacao)
Ts=1/fs;         % Periodo de Amostragem (s)
Ttotal=3;        % Tempo total de simulacao (s)
tempo=0:Ts:Ttotal;

% Parametros da planta - Filtro LCL (s) (parametros projetados de acordo 
% com a metodologia do artigo do Liserre

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
s = tf('s');
sysC3 = 1/(rc + rg + Lc*s + Lg*s + Cf*Lc*Lg*s^3 + Cf*Lc*rg*s^2 + ...
        Cf*Lg*rc*s^2 + Cf*rc*rg*s);
[numC3,denC3] = tfdata(sysC3,'v');

Lg_grid = 1e-3;       % indutancia da rede
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

% MR de 1 ordem 
% MR necessita ter ganho 0dB nas frequencias de interesse (nesse caso, 
% baixas frequencias) para nao dar um ganho no sinal de referencia
beta=8300;
alfa=-beta;
kMR=1;                        % ganho do modelo de referencia
numMR=kMR*[0 beta];
denMR=[1 -alfa];
sysMR=tf(numMR,denMR);

% Transformada Z, com ZOH, aplicada ao modelo de referencia Wm(s) 1 ordem
Gs4=tf(numMR,denMR);
Gz4=c2d(Gs4,Ts);
[numMRZ,denMRZ]=tfdata(Gz4,'v');
b0r=numMRZ(1);
b1r=numMRZ(2);
a0r=denMRZ(1);
a1r=denMRZ(2);

%% Parte 2 - Parametros PI-RMRAC

% Inicializacao dos vetores de entrada, saida e sinais internos (memoria)
y=zeros(1,length(tempo));           % Saida desejada
y1=zeros(1,length(tempo));          % Saida do modelo da planta - 1 ordem
y3=zeros(1,length(tempo));          % Saida do modelo da planta - 3 ordem
u=zeros(1,length(tempo));           % Lei de controle
ym=zeros(1,length(tempo));          % Saida do modelo de referencia
r=zeros(1,length(tempo));           % Sinal de referencia
e1=zeros(1,length(tempo));          % Erro de rastreamento
E1=zeros(1,length(tempo));          % Erro aumentado
m=zeros(1,length(tempo));           % Normalizador m
m2=4*ones(1,length(tempo));         % Normalizador m
mod_y=ones(1,length(tempo));        % Normalizador m
mod_u=ones(1,length(tempo));        % Normalizador m
vg=zeros(1,length(tempo));          % Tensao da rede (disturbio periodico)
d_cos=zeros(1,length(tempo));       % Decomposicao em fase e quadradatura
d_sin=zeros(1,length(tempo));       % do disturbio periodico (sin e cos)
norm_Theta=zeros(1,length(tempo));
sigma=zeros(1,length(tempo));
zetau=zeros*ones(1,length(tempo));   
zetauk1=zeros*ones(1,length(tempo));   
zetay=zeros*ones(1,length(tempo));   
zetae1=zeros*ones(1,length(tempo));   
zetaym=zeros*ones(1,length(tempo));   
zetac=zeros*ones(1,length(tempo));   
zetas=zeros*ones(1,length(tempo));   
Zeta=[zetau; zetauk1; zetay; zetae1; zetaym; zetac; zetas];

% Par�metros do controlador
Mo=20;           % Mo tem relacao com o valor da norma de Theta em RP
sigma_zero=0.1;  % Ts*lambda*sigma_k->1e-4*8000*0.125=0.1
delta0=0.7;
delta1=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Par�metros otimizados do controlador

% 5 agentes 100 itera��es
run1 = [11.837, 508.5441, -0.39627, -0.0064998, 0.3963, 0.69621, -0.078724, 0.019179, 0.11765];
run2 = [1, 500, -0.29339, 0.14184, 0.13454, 0.078213, 0.15965, 0.13897, 0.14753];
run3 = [3.5321, 598.912, -1.2971, -0.23385, 0.095052, 0.24889, -0.14707, 0.40109, 1.4538];
run4 = [1.2125, 715.3248, -0.37236, -0.048854, -0.038524, 0.25239, -1.2193, 0.039998, 0.47639];
run5 = [7.2651, 500, -1.2885, 0.7669, 0.10137, 0.70119, -0.43226, 0.16513, 0.14359];

% 10 agentes 100 itera��es
run6 = [4.7014, 660.8161, -0.20185, 0.10182, 0.20317, 0.17415, 0.18395, 0.14746, 0.052337, 25.4049];
run7 = [34.523, 578.5913, -1.7887, 0.052726, 0.034306, 0.92908, 0.5007, 0.1438, 1.6373, 0.11187];
run8 = [265.5619, 590.7392, -0.50666, -0.086886, 0.011716, 0.44834, -0.32646, 0.0056002, 0.49553, 0.074572];
run9 = [1.1803, 590.1398, 0.21153, 0.26905, -0.61153, 0.27962, 0.11286, -0.82041, -1.9131, 13335788.9462];
run10 = [1.0412, 520.934, -0.70352, 0.15328, 0.061639, 0.59095, 0.14452, 0.039015, 0.46886, 0.11779];

% 15 agentes 100 itera��es
run11 = [1.8541, 710.2661, -1.8246, 0.14274, -0.011546, 0.90747, 0.53606, 0.27579, 1.5939, 0.13065];
run12 = [4.0957, 892.9802, -1.749, 0.46496, -0.016117, 2, -0.095195, 0.052165, 1.1824, 0.068401];
run13 = [1.0024, 501.1768, -0.26747, 0.067298, -0.22757, 0.0026436, -0.28425, 0.00056861, 0.16587, 0.2847];
run14 = [4.0081, 578.7163, -0.60048, 0.076731, 0.171, 0.80034, -0.024933, -0.18412, -0.044732, 17.5527];
run15 = [218.7965, 558.9432, -0.8817, 0.40842, 0.0082336, 1.9225, -0.12118, 0.15772, 0.37372, 0.057107];

% 20 agentes 100 itera��es
run16 = [1.3099, 551.5609, -0.80244, -0.0019028, 0.46342, 0.82078, -0.55923, 0.13067, -0.42515, 38.046];
run17 = [3.3762, 571.683, -1.835, 0.88265, -0.22285, 0.19339, 0.20411, 0.23551, 0.86252, 0.3804];
run18 = [38.0116, 604.8815, -0.47799, -0.10762, 0.020136, 0.35052, -0.30885, -0.02343, 0.44618, 0.098003];
run19 = [41.4352, 662.8921, -0.92633, 0.10436, 0.059746, 1.5455, -0.89875, 0.025515, 0.23745, 0.10733];
run20 = [194.3537, 635.5666, -1.9901, 0.44327, 0.031362, 1.0557, 0.8953, 0.025369, 1.4343, 0.11453];

% 20 agentes 200 itera��es
run21 = [35.7515, 987.8531, -1.4072, 0.51663, -0.013731, 1.2659, -0.71557, 0.17614, 0.85765, 0.057765];
run22 = [2.1381, 500.301, -0.88779, 0.17626, 0.0042395, 0.75645, 0.051268, 0.040962, 0.62192, 0.061715];
run23 = [1.0288, 542.9777, -0.61523, 0.11242, -0.0062944, 1.1702, -0.75927, 0.13036, 0.46937, 0.054155];
run24 = [2.5701, 721.7404, -0.57441, -0.15811, 0.027216, 0.24692, 0.23961, -0.0037685, 0.57994, 0.15127];
run25 = [5.7335, 716.9691, -1.5274, 0.55911, 0.027067, 0.78496, -1.0734, 0.079632, 0.92519, 0.067329];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ganho = run18;

% Carrega par�metros e ganhos iniciais
Gamma = ganho(1);
kappa = ganho(2);
Theta1 = ganho(3) * ones(1, length(tempo));
Theta2 = ganho(4) * ones(1, length(tempo));
Theta3 = ganho(5) * ones(1, length(tempo));
Theta4 = ganho(6) * ones(1, length(tempo));
Theta5 = ganho(7) * ones(1, length(tempo));
Thetac = ganho(8) * ones(1, length(tempo));
Thetas = ganho(9) * ones(1, length(tempo));

Theta=[Theta1; Theta2; Theta3; Theta4; Theta5; Thetac; Thetas];

% Vari�veis para an�lise de desempenho
% dos regimes transit�rios
aux0 = 0;
aux1 = 0;
aux2 = 0;
aux3 = 0;
% dos overshoots/undershoots
ov0 = 0;
ov1 = 0;
ov2 = 0;
ov3 = 0;


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
    % e quadratura -> d_sin(k) e d_cos(k) ,  phi=pi/4;
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
    
    % Saida do Modelo de referencia - 1a. ordem
    ym(k)=-ym(k-1)*(a1r)+r(k-1)*(b1r);
    
    % Erro de rastreamento (e1)
    e1(k)=ym(k)-y(k);
    
    % Sigma modification
    norm_Theta(k)=norm(Theta(:,k));
    if (norm_Theta(k) < Mo)
        sigma(k)=0;
    elseif (norm_Theta(k) >=Mo && norm_Theta(k)<= 2*Mo)
        sigma(k)=sigma_zero*((norm_Theta(k)/Mo)-1);
    elseif (norm_Theta > 2*Mo)
        sigma(k)=sigma_zero;
    end
       
    % Vetor zeta (filtro auxiliar)
    zetau(k)  = -zetau(k-1)*(a1r)  + u(k-1)*(b1r);
    zetauk1(k)= -zetauk1(k-1)*(a1r)+ u(k-2)*(b1r);
    zetay(k)  = -zetay(k-1)*(a1r)  + y(k-1)*(b1r);
    zetae1(k) = -zetae1(k-1)*(a1r) + e1(k-2)*(b1r);
    zetaym(k) = -zetaym(k-1)*(a1r) + ym(k-1)*(b1r);
    zetac(k)  = -zetac(k-1)*(a1r)  + d_cos(k-1)*(b1r);
    zetas(k)  = -zetas(k-1)*(a1r)  + d_sin(k-1)*(b1r);
    
    Zeta(:,k) = [zetau(k);zetauk1(k);zetay(k);zetae1(k);zetaym(k);zetac(k);zetas(k)];
        
    % Erro aumentado
    E1(k)=y(k)+Theta(:,k)'*Zeta(:,k);

    % Lei de controle utilizada (descomentar a que for de interesse)
    u(k)=-(Theta2(k)*u(k-1)+Theta3(k)*y(k)+Theta4(k)*e1(k-1)+Theta5(k)*ym(k)+r(k))/Theta1(k)-(Thetac(k)*d_cos(k)+Thetas(k)*d_sin(k))/Theta1(k);
           
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

    m2(k) = (1-delta0*Ts)*m2(k-1) + delta1*Ts*( 1 + mod_y(k) + mod_u(k) );
    m(k)  = m2(k)*m2(k) + (Zeta(:,k)'*Gamma*Zeta(:,k));
    
    % Atualizacao dos parametros Theta (algoritmo gradiente modificado)
    Theta1(k+1)=Theta1(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*zetau(k)*E1(k))/([m(k)]);
    Theta2(k+1)=Theta2(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*zetauk1(k)*E1(k))/([m(k)]);
    Theta3(k+1)=Theta3(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*zetay(k)*E1(k))/([m(k)]);
    Theta4(k+1)=Theta4(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*zetae1(k)*E1(k))/([m(k)]);
    Theta5(k+1)=Theta5(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*zetaym(k)*E1(k))/([m(k)]);
    Thetac(k+1)=Thetac(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*zetac(k)*E1(k))/([m(k)]);
    Thetas(k+1)=Thetas(k)*(1-Ts*sigma(k)*Gamma)-(Ts*Gamma*kappa*zetas(k)*E1(k))/([m(k)]);
    Theta(:,k+1)=[Theta1(k+1);Theta2(k+1);Theta3(k+1);Theta4(k+1);Theta5(k+1);Thetac(k+1);Thetas(k+1)]; 
    
    %Transitorio Inicial
    if(tempo(k)>0 && tempo(k)<1)
        eval0 = abs(y(k));
        et0 = abs(e1(k));
        if(et0<=0.2 && aux0==0)
            t0 = tempo(k);
            aux0=1;
        end
        if(et0>0.2)
            aux0 = 0;
            %t0 = tempo(k);
        end
        if (eval0>ov0)
            ov0 = eval0;
        end 
    end

    
    %Transitorio de carga 1
    if(tempo(k)>=1 && tempo(k)<=1.75)
        eval1 = abs(y(k));
        et1 = abs(e1(k));
        if(et1<=0.4 && aux1==0)
            t1 = tempo(k)-1;
            aux1=1;
        end
        if(et1>0.4)
            aux1 = 0;
            %t1 = tempo(k)-1;
        end
        if (eval1>ov1)
            ov1 = eval1;
        end
    end
   
    
    %Transitorio de carga 2
    if(tempo(k)>1.75 && tempo(k)<=2.5)
        eval2 = abs(y(k));
        et2 = abs(e1(k));
        if(et2<=0.6 && aux2==0)
            t2 = tempo(k)-1.75;
            aux2=1;
        end
        if(et2>0.6)
            aux2 = 0;
            %t2 = tempo(k)-1.75;
        end
        if (eval2>ov2)
            ov2 = eval2;
        end
    end

  
    %Transitorio de variacao parametrica
    if(tempo(k)>2.5)
        eval3 = abs(y(k));
        et3 = abs(e1(k));
        if(et3<=0.6 && aux3==0)
            t3 = tempo(k)-2.5;
            aux3=1;
        end
        if(et3>0.6)
            aux3 = 0;
            %t3 = tempo(k)-2.5;
        end
        if (eval3>ov3)
            ov3 = eval3;
        end
    end

      
    
end

% Calcula m�tricas do erro
MAE = mad(e1);
MSE = immse(r,y);
RMSE = sqrt(mean((r - y).^2));


% Overshoots/Undershoots
    % inicial
    if (ov0<10)
        ov0 = 0;
    else
        ov0 = ov0 - 10;
    end
    
    % carga 1
    if (ov1<20)
        ov1 = 0;
    else
        ov1 = ov1 - 20;
    end
    
    % carga 2
    if (ov2<30)
        ov2 = 0;
    else
        ov2 = ov2 - 30;
    end
    
    % varia��o param�trica
    if (ov3<30)
        ov3 = 0;
    else
        ov3 = ov3 - 30;
    end

% Documenta �ndices de performance
logFilename = 'Outputs/performance.txt';
diary(logFilename); 

display(['MAE: ', num2str(MAE)]);
display(['MSE: ', num2str(MSE)]);
display(['RMSE: ', num2str(RMSE)]);
display(['T0: ', num2str(t0)]);
display(['T1: ', num2str(t1)]);
display(['T2: ', num2str(t2)]);
display(['T3: ', num2str(t3)]);
display(['Overshoot 0: ', num2str(ov0)]);
display(['Overshoot 1: ', num2str(ov1)]);
display(['Overshoot 2: ', num2str(ov2)]);
display(['Overshoot 3: ', num2str(ov3)]);
display(' ');

diary off;


%% Parte 4 - Graficos 

fonte = 21;

a=figure;
% plot(tempo,r,'k','LineWidth',2);
plot(tempo,ym,':k','LineWidth',4);
hold on
% plot(t,y,'color',[0,0,0]+0.7,'LineWidth',2);
plot(tempo,y,'b','LineWidth',2); %4
grid on;
h2 = legend('$y_m$','$y$',([345, 100, 0, 0]));
% h2 = legend('$y_{m}$','$y$',([345, 100, 0, 0]));
set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
set(gcf,'Units','centimeters','PaperSize',[13 7]);   
% Recortar a figura da p�gina
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
plot(tempo,e1,'b','LineWidth',3);
hold
% plot(tempo,E1,'color',[0,0,0]+0.7,'LineWidth',2);
% plot(tempo,E1,':k','LineWidth',4);
grid on;
% h2 = legend('$e_{1}$','$\epsilon$',([345, 100, 0, 0]));
h2 = legend('$e_{1}$',([345, 100, 0, 0]));
set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
set(gcf,'Units','centimeters','PaperSize',[13 7]);                                     % Recortar a figura da p�gina
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
set(gcf,'Units','centimeters','PaperSize',[13 7]);                                     % Recortar a figura da p�gina
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
plot(tempo,Theta(7,:),'--k','LineWidth',2);
grid on;
h2 = legend('$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$','$\theta_{5}$','$\theta_{c}$','$\theta_{s}$',([345, 100, 0, 0]));
set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
set(gcf,'Units','centimeters','PaperSize',[14 7]);                                     % Recortar a figura da p�gina
set(gca,'fontsize',fonte,'units','norm');
ylabel('Gains','fontsize',fonte);
xlabel('Time (s)','fontsize',fonte);
% ylim([-4 8])