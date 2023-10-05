function res = PI_RMRAC_for_optimization(x)

    %% Parte 1 - Definicoes do sistema

    fs=5.04e3;              % Frequencia de discretizacao ( = Passo da simulacao)
    f=60;                   % Frequencia do sinal de interesse
    Ts=1/fs;                % Periodo de Amostragem (s)
    Ttotal=3;               % Tempo total de simulacao (s)
    tempo=0:Ts:Ttotal;

    Vlink = 400;   % tensao de barramento (monofasico 127V rede)
    Cf = 62e-6;    % capacitor filtro LCL

    rc = 50e-3;    % resistencia referente ao primeiro braco do LCL
    rg = 50e-3;    % resistencia referente segundo braco do LCL

    Rd = 1;        % resistencia de amortecimento passivo (serie capacitor LCL)
    Ro = 0;        % curto-circuito terminais fonte (sem considerar fonte Vg)

    Lc = 1e-3;     % Indutancia referente ao primeiro braco do LCL
    Lg = 0.3e-3;   % Indutancia referente ao segundo braco do LCL

    % Modelo de terceira ordem - para controle da corrente do lado da rede (ig)
    s = tf('s');
    sysC3 = 1/(rc + rg + Lc*s + Lg*s + Cf*Lc*Lg*s^3 + Cf*Lc*rg*s^2 + ...
            Cf*Lg*rc*s^2 + Cf*rc*rg*s);
    [numC3,denC3] = tfdata(sysC3,'v');

    Lg_grid = 1e-3;       % indutancia da rede (variação paramétrica)
    
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
     

    %% PI-RMRAC
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
    
    delta0=0.7;
    delta1=1;  
    Mo = 20;        % Mo tem relacao com o valor da norma de Theta em RP
    sigma_zero=0.1; % Ts*lambda*sigma_k->1e-4*8000*0.125=0.1
   
    % Parâmetros otimizados do controlador
    
    Gamma=x(1); 
    kappa=x(2);
    
    % Valor Inicial dos Ganhos de adaptacao

    Theta1=x(3)*ones(1,length(tempo));      
    Theta2=x(4)*ones(1,length(tempo)); 
    Theta3=x(5)*ones(1,length(tempo)); 
    Theta4=x(6)*ones(1,length(tempo)); 
    Theta5=x(7)*ones(1,length(tempo)); 
    Thetac=x(8)*ones(1,length(tempo));    
    Thetas=x(9)*ones(1,length(tempo)); 
    
    Theta=[Theta1; Theta2; Theta3; Theta4; Theta5; Thetac; Thetas];
    
    % Simulação
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
        % e quadratura -> d_sin(k) e d_cos(k) , phi=pi/4;
        phi=0;
        
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
        ym(k)= -ym(k-1)*(a1r)+r(k-1)*(b1r);

        % Erro de rastreamento (e1)
        e1(k)= ym(k)-y(k);

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

        Zeta(:,k)=[zetau(k);zetauk1(k);zetay(k);zetae1(k);zetaym(k);zetac(k);zetas(k)];

        % Erro aumentado
        E1(k)=y(k)+Theta(:,k)'*Zeta(:,k);

        % Lei de controle utilizada (descomentar a que for de interesse)
        u(k)=-(Theta2(k)*u(k-1)+Theta3(k)*y(k-1)+Theta4(k)*e1(k-1)+Theta5(k)*ym(k)+r(k))/Theta1(k)-(Thetac(k)*d_cos(k)+Thetas(k)*d_sin(k))/Theta1(k);

        % Modulos e funcao normalizacao
        if u(k-1) < 0
            mod_u(k)=-u(k);
        else
            mod_u(k)=u(k);        
        end

        if y(k-1) < 0
            mod_y(k)=-y(k);
        else
            mod_y(k)=y(k);    
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
            
    end

 
    %% Penalizações da função custo
    o = mean(abs(e1));

    % testa valores globais    
    if max(y) > 40 || max(e1) > 40 || max(u) > 400 || o==NaN
        o=o*10;            
    end    
    if min(y) < -40 || min(e1) < -40 || min(u) < -400
        o=o*10;            
    end
     
    % restricao prova de estabilidade do PI
    if 2*Ts*Gamma*kappa*sigma(k)*sigma(k)*(norm(Zeta(:,k))*Zeta(:,k))/m2(k) > 0.99
       o=o*20;
    end
    
    
end