clear all 
close all
clc

Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper
SearchAgents_no=30; % Number of search agents
Max_iteration=500; % Maximum numbef of iterations

num = 1; %Numero de vezes que o codigo roda
logFilename = 'ganhos.txt';
imgCounter = 1;
diary(logFilename); 

for run = 1:num

tic
% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
[Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
toc 


    display(['Run: ', num2str(run)]);
    display(['Ganhos: ', num2str(Best_pos)]);
    display(['Fitness : ', num2str(Best_score)]);

    
    fonte = 21;
    figure
    hold on
    semilogy(GWO_cg_curve,'Color','b','LineWidth',3);
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
    filename = "fitness_"+int2str(imgCounter)+".pdf";
    exportgraphics(gcf, filename, 'ContentType', 'vector', 'BackgroundColor', 'none');
    imgCounter = imgCounter+1;
    
end
diary off;
