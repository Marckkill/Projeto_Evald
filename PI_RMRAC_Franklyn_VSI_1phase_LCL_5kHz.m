%========================================================================
% Universidade Federal de Santa Maria
% Guilherme Hollweg
% 
% Algoritmo PI-RMRAC aplicado na planta do VSI com LCL monofasico 
%
% Data: 31/03/2022
%
% Ultima Modificacao: Paulo Evald(31/03/2022)
%========================================================================

close all
clear all
clc

format long g

%% Parte 0 - Definicoes para os Diagramas de Bode

% paramentros de configuracao dos diagramas de Bode
opts=bodeoptions('cstprefs');
opts.Title.Interpreter='latex';
opts.Title.FontSize=12;
opts.XLabel.Interpreter='latex';
opts.XLabel.String='Frequ\^encia';
opts.XLabel.FontSize=12;
opts.YLabel.Interpreter='latex';
opts.YLabel.FontSize=12;
opts.FreqUnits='Hz';
opts.Grid='on';

%% Parte 1 - Parametros de Simulacao e Definicao Plantas
fs=5.04e3;       % Frequencia de discretizacao ( = Passo da simulacao)
Ts=1/fs;         % Periodo de Amostragem (s)
Ttotal=3;        % Tempo total de simulacao (s)
tempo=0:Ts:Ttotal;

% Parametros da planta - Filtro LCL (s) (parametros projetados de acordo 
% com 'Liserre', e simulados no PSIM para VSI com Controle Classico)
%
% Necessario que os resistores 'rc' e 'rg' nao tenham valor nulo, pois o 
% pico de ressonancia (sem amortecimento) fica muito complicado de ser
% controlado, e o controle facilmente tende a instabilidade

Vlink = 500;   % tensao de barramento (monofasico 127V rede)
Cf = 62e-6;    % capacitor filtro LCL

%rc = 0;
%rg = 0;
rc = 50e-3;    % resistencia referente ao primeiro braco do LCL
rg = 50e-3;    % resistencia referente segundo braco do LCL

%Rd = 1;        % resistencia de amortecimento passivo (serie capacitor LCL)
Rd = 0;
Ro = 0;        % curto-circuito terminais fonte (sem considerar fonte Vg)

Lc = 1e-3;     % Indutancia referente ao primeiro braco do LCL
Lg = 0.3e-3;   % Indutancia referente ao segundo braco do LCL

% Modelo de terceira ordem - para controle da corrente do lado da rede (ig)
% - conforme a Tese - Anexo C - Modelagem da planta
s = tf('s')
sysC3 = 1/(rc + rg + Lc*s + Lg*s + Cf*Lc*Lg*s^3 + Cf*Lc*rg*s^2 + ...
        Cf*Lg*rc*s^2 + Cf*rc*rg*s)
[numC3,denC3] = tfdata(sysC3,'v');

Lg_grid = 1e-3;       % indutancia da rede
% Modelo da planta completa considerando a variacao parametrica
sysC3_VP = 1/(rc + rg + Lc*s + (Lg+Lg_grid)*s + Cf*Lc*(Lg+Lg_grid)*s^3 + Cf*Lc*rg*s^2 + ...
        Cf*(Lg+Lg_grid)*rc*s^2 + Cf*rc*rg*s)
[numC3_VP,denC3_VP] = tfdata(sysC3_VP,'v');    

% Modelo considerando Rd
sysC3_Rd = ((Cf * Rd * s + 1) ) / (rc + rg + Ro + Lc * s ...
        + Lg * s + Cf * Lc * Lg * s^3 + Cf * Lc * rg * s^2 + Cf * Lg ...
         * rc * s^2 + Cf * Lc * Rd * s^2 + Cf * Lg * Rd * s^2 + ...
         Cf * Lc * Ro * s^2 + Cf * rc * rg * s + Cf * rc * rg * s ...
         + Cf * rg * Rd * s + Cf * rc * Ro * s + Cf * Rd * Ro * s)

% Gp(s) - Parte Modelada da Planta - 1 ordem
% negligenciando os polos complexos conjugados
kpc = 1 / (Lc + Lg);
numC = kpc;
denC = [1 (rg + rc) / (Lc + Lg) ];
sysC = tf( numC , denC)

% Polo real do modelo de primeira ordem
p4 = abs( roots (denC) )

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
% Mesmo MR da Tese, consultar Anexo D
beta=8300;
alfa=-beta;
% beta=30000;
% alfa=-beta;
kMR=1;                        % ganho do modelo de referencia
numMR=kMR*[0 beta];
denMR=[1 -alfa];
sysMR=tf(numMR,denMR)

% Transformada Z, com ZOH, aplicada ao modelo de referencia Wm(s) 1 ordem
Gs4=tf(numMR,denMR);
Gz4=c2d(Gs4,Ts);
[numMRZ,denMRZ]=tfdata(Gz4,'v');
b0r=numMRZ(1)
b1r=numMRZ(2)
a0r=denMRZ(1);
a1r=denMRZ(2)

delta0=0.7;
delta1=1;

% figure(1)
% b2=bodeplot(sysC3, sysC, opts);
% title('Diagrama de Bode Planta Real e Aproximada')
% l2=legend('$G\_{(ig,d)}(s) Real$', '$G\_{(ig,d)}(s) Aprox.$'); 
% set(l2, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% axes=findobj('type','axes');
% b2_magnitude=get(axes(2),'YLabel');
% b2_phase=get(axes(1),'YLabel');
% set(b2_magnitude,'String','Magnitude (dB)');
% set(b2_phase,'String','Fase (graus)');
% 
% figure(2)
% b3=bodeplot(sysC, sysMR, opts);
% title('Diagrama de Bode Planta Aproximada e Modelo Refer\^encia')
% l3=legend('$Gp(s) Aprox.$', '$Wm(s)$'); 
% set(l3, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% axes=findobj('type','axes');
% b3_magnitude=get(axes(2),'YLabel');
% b3_phase=get(axes(1),'YLabel');
% set(b3_magnitude,'String','Magnitude (dB)');
% set(b3_phase,'String','Fase (graus)');

%% Parte 2 - Parametros PI-RMRAC

% Parametros da lei de adaptacao parametrica

% Ganhos obtidos com GA
% Gamma=801;                             % normalizador
% kappa=997.601423622602;                % ganho de adaptacao thetas

% Ganhos que eu sintonizei
% Gamma=1000;                 % normalizador
% kappa=3;                 % ganho de adaptacao thetas

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
% sigma_zero=0.5;                   % Ts*lambda*sigma_k->1e-4*8000*0.125=0.1
% Mo=0.2/2;                         % Mo tem relacao com o valor da norma
                                  % de Theta em RP, podendo ser metade
                                  % do valor, tendo em vista que a sigma
                                  % modification usa Mo e 2*Mo

sigma_zero=0.1;
Mo=5;

% Valor Inicial dos Ganhos de adaptacao

% Aleatório
% Gamma = 1;
% kappa = 1000;
% Theta1 = -5*ones(1,length(tempo));
% Theta2 = 1*ones(1,length(tempo));
% Theta3 = -2*ones(1,length(tempo));
% Theta4 = 1*ones(1,length(tempo));
% Theta5 = 1*ones(1,length(tempo));
% Thetac = 0.1*ones(1,length(tempo));
% Thetas = 0.1*ones(1,length(tempo));
% Theta=[Theta1; Theta2; Theta3; Theta4; Theta5; Thetac; Thetas];

% com k=4000
% Theta1 = -3.69596652301644*ones(1,length(tempo));
% Theta2 = 2.21122338417208*ones(1,length(tempo));
% Theta3 = -2.60095474560142*ones(1,length(tempo));
% Theta4 = 1.66506624984036*ones(1,length(tempo));
% Theta5 = 1.06642094694897*ones(1,length(tempo));
% Thetac = 0.309372903825546*ones(1,length(tempo));
% Thetas = 1.51137460586707*ones(1,length(tempo));
% Theta=[Theta1; Theta2; Theta3; Theta4; Theta5; Thetac; Thetas];

run1 = [213.0443, 239.375, -0.6755966, 0.1534227, 0.02251536, 0.6236674, -1.119297, 0.03257845, 0.5142193];
run2 = [101.5284, 759.2285, -0.5872466, 0.2532434, -0.02019455, 0.8453407, -0.7960242, -0.02258432, 0.2275902];
run3 = [205.93, 207.2881, -0.772152, 0.5102774, 0.01055886, 0.6724002, -0.4263753, 0.04409399, 0.2033881];
run4 = [196.5503, 441.2579, -0.4855506, 0.1405096, 0.01270843, 0.8899538, 0.2440502, -0.02063314, 0.2345386];
run5 = [115.3526, 115.462, -0.5046487, 0.2898377, -0.08364972, 0.129221, -0.1305447, 0.02118995, 0.1667776];
run6 = [53.38864, 161.0895, -0.5286447, 0.1045305, 0.0210228, 1.146222, -1.050839, 0.03277853, 0.3442126];
run7 = [8.0845, 8.1467, -0.84418, 0.43163, -0.97496, -0.36209, -0.39496, 0.058619, 0.43048]; % ERRO
run8 = [124.9423, 127.9772, -0.4286573, 0.08442583, -0.01317953, 0.231005, 0.6017969, 0.08958855, 0.2555118];
run9 = [217.8912, 226.6223, -0.3279347, 0.1452973, -0.07297243, 0.2277312, -0.3341671, -0.009822892, 0.1420925];
run10 = [172.956, 245.6202, -5.033108, 3.62641, 0.04974392, 4.595189, 3.276951, 0.3981488, 1.12355];

ganho = run10;

% GA 
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

zetau=zeros*ones(1,length(tempo));   
zetauk1=zeros*ones(1,length(tempo));   
zetay=zeros*ones(1,length(tempo));   
zetae1=zeros*ones(1,length(tempo));   
zetaym=zeros*ones(1,length(tempo));   
zetac=zeros*ones(1,length(tempo));   
zetas=zeros*ones(1,length(tempo));   
Zeta=[zetau; zetauk1; zetay; zetae1; zetaym; zetac; zetas];
%% Parte 3 - Simulacao do Sistema em Malha Fechada (PI-RMRAC)
aux0 = 0;
aux1 = 0;
aux2 = 0;
aux3 = 0;

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
    
    % Atualizacao dos parametros Theta (algoritmo gradiente)
    % solucao completa (homogenea + particular)
    % usar erro aumentado 
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
        et0 = abs(e1(k));
        if(et0<=0.2 && aux0==0)
            t0 = tempo(k);
            aux0=1;
        end
        if(et0>0.2)
            aux0 = 0;
            %t0 = tempo(k);
        end
    end

    
    %Transitorio de carga 1
    if(tempo(k)>=1 && tempo(k)<=1.75)
        et1 = abs(e1(k));
        if(et1<=0.4 && aux1==0)
            t1 = tempo(k)-1;
            aux1=1;
        end
        if(et1>0.4)
            aux1 = 0;
            %t1 = tempo(k)-1;
        end
    end
    


    
    %Transitorio de carga 2
    if(tempo(k)>1.75 && tempo(k)<=2.5)
        et2 = abs(e1(k));
        if(et2<=0.6 && aux2==0)
            t2 = tempo(k)-1.75;
            aux2=1;
        end
        if(et2>0.6)
            aux2 = 0;
            %t2 = tempo(k)-1.75;
        end
    end


    
    
    
    %Transitorio de variacao parametrica
    if(tempo(k)>2.5)
        et3 = abs(e1(k));
        if(et3<=0.6 && aux3==0)
            t3 = tempo(k)-2.5;
            aux3=1;
        end
        if(et3>0.6)
            aux3 = 0;
            %t3 = tempo(k)-2.5;
        end
    end
end

MAE = mad(e1)
MSE = immse(r,y)
RMSE = sqrt(mean((r - y).^2))


display(t0);
display(t1);
display(t2);
display(t3);

%% Parte 4 - Graficos 

% figure(3)
% plot(tempo,r);
% grid on
% title('Refer\^encia', 'Interpreter', 'latex', 'FontSize', 16);
% l4=legend('r');
% set(l4, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% xlabel('Tempo (s)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Amplitude Refer\^encia r', 'Interpreter', 'latex', 'FontSize', 14);
% 
% figure(4)
% plot(tempo,e1,'r');
% hold on
% plot(tempo,E1,'b');
% grid on
% title('Erros Rastreamento x Aumentado', 'Interpreter', ...
%       'latex', 'FontSize', 16);
% l6=legend('e1 (rastreamento)', 'E1 (aumentado)');
% set(l6, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% xlabel('Tempo (s)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Amplitude dos Erros','Interpreter', ...
%        'latex', 'FontSize', 14);
%    
% figure(5)
% plot(tempo,u);
% grid on
% title('Acao de Controle Escolhida', 'Interpreter', ... 
%       'latex', 'FontSize', 16);
% l8=legend('u');
% set(l8, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% xlabel('Tempo (s)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Amplitude da Lei de Controle (MRAC+SM)','Interpreter', ...
%        'latex', 'FontSize', 14);
% 
% figure(6)
% plot(tempo,Theta1(1,1:end-1),'r');
% hold on
% plot(tempo,Theta2(1,1:end-1),'b');
% hold on
% plot(tempo,Theta3(1,1:end-1),'g');
% hold on
% plot(tempo,Theta4(1,1:end-1),'y');
% hold on
% plot(tempo,Theta5(1,1:end-1),'c');
% hold on
% plot(tempo,Thetac(1,1:end-1),'k');
% hold on
% plot(tempo,Thetas(1,1:end-1),'r');
% grid on
% title('Ganhos Algoritmo Adaptacao (Gradiente)', ...
%       'Interpreter', 'latex', 'FontSize', 16);
% l9=legend('Theta1', 'Theta2', 'Theta3', 'Theta4', 'Theta5', 'ThetaC', 'ThetaS');
% set(l9, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% xlabel('Tempo (s)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Amplitude de Ganhos (Theta)','Interpreter', ...
%        'latex', 'FontSize', 14);
% 
% figure(7)
% plot(tempo,m);
% ylabel('m');
% grid on;
% title('Normalizador', 'Interpreter', 'latex', 'FontSize', 16);
% l10=legend('m');
% set(l10, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% xlabel('Tempo (s)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Amplitude do Normalizador','Interpreter', ...
%        'latex', 'FontSize', 14);
% 
% figure(8)
% plot(tempo,ym,'r');
% hold on
% plot(tempo,y,'b');
% grid on
% title('Saidas do sistema', 'Interpreter', 'latex', 'FontSize', 16);
% l11=legend('ym', 'y');
% set(l11, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% xlabel('Tempo (s)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Amplitude de saida', 'Interpreter', 'latex', 'FontSize', 14);
% 
% figure(9)
% plot(tempo,norm_Theta,'r');
% grid on
% title('Norma de Theta', 'Interpreter', 'latex', ...
%       'FontSize', 16);
% l12=legend('Norma de Theta');
% set(l12, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% xlabel('Tempo (s)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Norma de Theta','Interpreter', 'latex', ...
%        'FontSize', 14);
%    
% figure(10)
% plot(tempo,sigma,'r');
% grid on
% title('Sigma Modification', 'Interpreter', 'latex', ...
%       'FontSize', 16);
% l13=legend('sigma');
% set(l13, 'Interpreter', 'latex', 'FontSize', 10, 'Orientation', ...
%     'vertical', 'location', 'northeast');
% xlabel('Tempo (s)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Sigmas','Interpreter', 'latex', ...
%        'FontSize', 14);
% 


% PI-RMRAC plots para artigo
fonte = 21;

a=figure;
% plot(tempo,r,'k','LineWidth',2);
plot(tempo,ym,'k','LineWidth',2);
hold on
% plot(t,y,'color',[0,0,0]+0.7,'LineWidth',2);
plot(tempo,y,':b','LineWidth',4); %4
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
plot(tempo,E1,':k','LineWidth',4);
grid on;
h2 = legend('$e_{1}$','$\epsilon$',([345, 100, 0, 0]));
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

% t(k+1) = t(k);
% figure
% set(gcf,'renderer','painters')
% plot(t,theta(1,:),':k','LineWidth',4);
% hold on
% plot(t,theta(2,:),'color',[0,0,0]+0.6,'LineWidth',3);
% plot(t,theta(3,:),'--k','LineWidth',2);
% plot(t,theta(4,:),'-.k','LineWidth',2);
% plot(t,theta(5,:),'-.','color',[0,0,0]+0.7,'LineWidth',2);
% grid on;
% h2 = legend('$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$','$\theta_{5}$',([345, 100, 0, 0]));
% set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
% set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
% set(gcf,'Units','centimeters','PaperSize',[14 7]);                                     % Recortar a figura da p�gina
% set(gca,'fontsize',fonte,'units','norm');
% ylabel('Gains','fontsize',fonte);
% xlabel('Time (s)','fontsize',fonte);

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