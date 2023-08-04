clear all %#ok<CLALL>
close all
clc

N=20; %Numero de falcoes
Function_name='F1'; %Sistema  
T=200; %Numero de iteracoes

num = 10; %Numero de vezes que o codigo roda
logFilename = 'console_output_log.txt';
diary(logFilename);

for run = 1:num
    [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
    [Rabbit_Energy,Rabbit_Location,CNVG]=HHO(N,T,lb,ub,dim,fobj);
    display(['Loop: ', num2str(run)]);
    display(['Numero de falcoes: ', num2str(N)]);
    display(['Numero de iteracoes: ', num2str(T)]);
    display(['Valores otimizados: ', num2str(Rabbit_Location)]);
    display(['Melhor fitness: ', num2str(Rabbit_Energy)]);
    display(' ');
end

diary off;