clear all
close all
clc

N=50; % Number of search agents
T=500; % Maximum number of iterations
Function_name='F1'; % Name of the test function

outputFolder = 'Outputs_Final'; %Nome da pasta de saida
imgCounter = 1; %Contador de figuras
num = 5; %Numero de vezes que o codigo roda
logFilename = 'ganhos.txt'; %Nome do log
logFullfile = fullfile(outputFolder,logFilename);

for run = 1:num
    
    
    % Load details of the selected benchmark function
    [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
    [Rabbit_Energy,Rabbit_Location,CNVG]=HHO(N,T,lb,ub,dim,fobj);
    
    diary(logFullfile); %Inicia o log
    % Display information
    display(['Run: ', num2str(run)]);
    display(['Ganhos: ', num2str(Rabbit_Location)]);
    display(['Mean absolute error: ', num2str(Rabbit_Energy)]);
    fprintf('\n');
    diary off;
    
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
    
    
    % Salva a figura
    set(gcf, 'Position', get(0, 'ScreenSize'));
    filename = fullfile(outputFolder, ['fitness_', num2str(imgCounter), '.pdf']);
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    imgCounter = imgCounter + 1;
    
    
    
end
