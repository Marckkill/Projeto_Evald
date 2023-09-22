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
run11 = [419.6988, 416.0337, -0.297972, 0.1420171, 0.1103898, 0.09699574, 0.1528626, 0.1596322, 0.1450009];
run12 = [167.1965, 75.1796, -1.53409, 1.267311, 0.05832269, 0.2874038, -0.4156018, 0.03451441, 0.05855093];
run13 = [250.4957, 229.2439, -0.3013227, -0.0002652771, 0.09490825, 0.1763071, 0.05000762, -0.08538051, 9.704084e-05];
run14 = [124.3937, 124.3946, -0.1252888, 0.04785739, -0.02437564, -0.02211184, -0.0324554, -0.03978009, -0.05963512];
run15 = [70.68685, 194.6275, -0.1830204, 0.03927152, 0.02798445, 0.3346562, -0.05516622, -0.02336735, -0.02146523];
run16 = [220.1712, 535.8073, -0.4556295, 0.1605824, 0.1636634, 0.2456095, 1.424088, 0.282578, 0.200289];
run17 = [423.4247, 630.5384, -0.6307056, 0.1032169, 0.2878185, 0.6591699, 0.4457106, 0.5667509, 0.09940332];
run18 = [85.2097, 64.0382, -0.172098, 0.0475738, -0.200376, -0.0543612, -0.068621, 0.00909796, 0.0813793]; %ERRO
run19 = [305.9976, 319.9186, -0.2395501, 0.07405938, 0.01475971, 0.1540682, 0.4546484, 0.1181367, 0.1464548];
run20 = [74.62692, 113.3649, -0.4039666, 0.0470856, 0.02169757, 0.3729628, 0.02453606, 0.0356503, 0.2830511];
run21 = [254.056, 589.8777, -0.956227, 0.1470321, 0.006317114, 0.4934036, 0.1857939, 0.1589699, 0.7548527];
run22 = [826.5501, 198.2087, -0.2965418, 0.005248478, 0.09897313, 0.1350394, 0.1272492, -0.05301344, -0.009200613];
run23 = [796.0095, 729.3894, -1.082994, 0.4232208, 0.03025466, 0.5974452, 0.5263238, -0.01671684, 0.4731746];
run24 = [106.7286, 151.8981, -0.1440115, 0.05045504, -0.1148311, 0.04216321, 0.04095441, 0.007005651, 0.04419522];
run25 = [318.5635, 312.9175, -0.260408, 0.2132526, 0.1120209, 0.2191224, 0.2054765, 0.1059433, -0.01764737];
run26 = [209.5551, 209.4442, -0.1316146, 0.04589244, 0.08714409, 0.08780425, -0.9561599, 0.1004829, 0.08395257];
run27 = [878.7223, 618.6148, -0.670339, 0.0183705, 0.04130391, 0.8936983, 0.6312078, -0.008736371, 0.1629408];
run28 = [25.19652, 196.3119, -0.192172, 0.05705789, 0.006116985, -0.006573285, -0.05722085, -0.09198096, 0.002113664];
run29 = [1, 1, -1.632, 0.65932, -1.5236, 3.4339, 3.4483, -1.7401, 1.8142];
run30 = [366.2374, 438.35, -0.91075, 0.3693762, 0.006617532, 0.644127, -0.1776385, 0.1089609, 0.5119386];

% kappa ub = 500
run31 = [1, 500, 4.204698, 5.157824, -2.213767, -4.017727, 5.520255, -2.60423, -4.986212];
run32 = [1.111651, 645.0597, -0.5676039, 0.2321534, 0.1780938, 0.7219489, -0.175477, -0.05033073, -0.1619689];
run33 = [1.003738, 500, -0.2065758, 0.1567696, 0.2201972, 0.2696581, -0.1168999, -0.12811, -0.118407];
run34 = [1.009529, 504.8071, -5.815949, -1.953983, 7.985341, -3.792057, -4.57856, -3.884587, 1.047199];
run35 = [1.012757, 500, -12.62469, -0.6160999, -2.274302, -13.01796, -2.232711, -16.78853, -15.70147];
run36 = [1.516513, 753.7504, 4.81704, 4.586132, 4.337695, 4.577723, 4.557078, 4.747633, 4.657636];
run37 = [1.690616, 550.238, -4.762017, 2.033934, 0.5728005, 1.106657, 2.35903, 0.8199304, 1.18093];
run38 = [1.999145, 575.3971, -1.21653, -0.4216469, -0.01010372, 0.1874922, -0.01872024, 0.3381104, 2.084742];
run39 = [2.936327, 872.1612, -2.333776, 0.8039738, 0.02770514, 1.785212, 0.5163204, 0.3773891, 1.47251];
run40 = [1.000578, 500.2889, -13.43323, 0.3468666, -0.494856, -11.76697, -13.86063, -2.030159, -13.92134];
run41 = [3.275421, 504.3389, -1.320305, 0.4435987, -0.1122919, 1.207823, 0.09387756, 0.1491838, 0.8860272];
run42 = [1.667531, 520.7365, -0.5737836, 0.1900324, 0.2199819, 0.405078, -0.5966206, -0.1237136, -0.1220482];
run43 = [3.548399, 856.7041, -0.5743915, 0.2338526, 0.06820062, 0.4904971, 0.1056865, -0.1139974, 0.06484549];
run44 = [7.372292, 513.8119, -1.584642, -1.470564, 0.1102533, 0.5291499, -2.255399, 0.1700786, 3.136364];
run45 = [1.53783, 681.8151, -0.7758083, 0.1892139, 0.03443082, 1.51919, -0.4295469, -0.2732915, 0.4375035];
run46 = [1.485978, 605.4617, -0.611414, 0.3056595, 0.01457811, 0.5406216, -0.01046334, -0.03677821, 0.1698201];
run47 = [2.520014, 506.935, -0.6624829, 0.3146231, 0.01466696, 0.4776982, -0.007233049, 0.1355761, 0.2909992];
run48 = [21.82886, 502.3119, -1.404129, 0.2884015, 0.1025483, 1.270097, 0.2512866, 0.0606212, 0.1794631];
run49 = [17.09666, 524.4375, -4.294781, 3.061051, 0.007783136, 2.345235, -1.495609, 0.2651109, 1.186419];
run50 = [1.007286, 503.6428, -0.2953276, 0.04297083, 0.07180788, 0.2148896, -0.04200504, -0.1235946, 0.1082169];





ganho = run47;

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
ov0 = 0;
ov1 = 0;
ov2 = 0;
ov3 = 0;

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

if (ov0<10)
    ov0 = 0;
else
    ov0 = ov0 - 10;
end

if (ov1<20)
    ov1 = 0;
else
    ov1 = ov1 - 20;
end

if (ov2<30)
    ov2 = 0;
else
    ov2 = ov2 - 30;
end

if (ov3<30)
    ov3 = 0;
else
    ov3 = ov3 - 30;
end


MAE = mad(e1);
MSE = immse(r,y);
RMSE = sqrt(mean((r - y).^2));

logFilename = 'erros.txt';
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

hold off