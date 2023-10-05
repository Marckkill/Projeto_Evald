clear all 
close all
clc

N=30; % Number of search agents
T=500; % Maximum number of iterations
Function_name='F1'; % Name of the test function 

num = 1; %Numero de vezes que o codigo roda
logFilename = 'ganhos.txt';
imgCounter = 1;
diary(logFilename); 
for run = 1:num

    % Load details of the selected benchmark function
    [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
    [Rabbit_Energy,Rabbit_Location,CNVG]=HHO(N,T,lb,ub,dim,fobj);

    % Display information
    display(['Best gains are: ', num2str(Rabbit_Location)]);
    display(['Mean absolute error: ', num2str(Rabbit_Energy)]);

    % Generate figure of fitness convergence
    fonte = 21;
    figure
    hold on
    semilogy(CNVG,'Color','b','LineWidth',3);
    h2 = legend('Fitness',([345, 100, 0, 0]));
    set(h2, 'interpreter','latex','fontsize',fonte,'units','norm','Location','NorthEast'); % Legenda
    set(gcf,'Units','centimeters','Position',[10,7,8.8,5.3],'color','white');              % Background
    set(gcf,'Units','centimeters','PaperSize',[13 7]);                                     % Recortar a figura da página
    set(gca,'fontsize',fonte,'units','norm');
    xlabel('Iteration','fontsize',fonte);
    ylabel('Mean absolute error','fontsize',fonte);
    axis tight
    grid off
    box on

    saveas(gcf,"fitness_"+int2str(imgCounter)+".pdf")
    imgCounter = imgCounter+1;

end
diary off;
