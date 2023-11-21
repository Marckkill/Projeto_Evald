%_________________________________________________________________________%
%  Black Widow Optimization Algorithm (source code)                       %
%                                                                         %
%  Developed in MATLAB R2018a(9.4)                                        %
%                                                                         %
%  Author:  Dr. Hernan Peraza-Vazquez                                     %
%                                                                         %
%  Programmer: Hernan Peraza-Vazquez                                      %
%                                                                         %
%         e-Mails: hperaza@ipn.mx   or   hperaza@ieee.org                 %
%                                                                         %
%                                                                         %
% Paper:  A Novel Bio-Inspired Algorithm Applied to Selective Harmonic    %
%         Elimination in a Three-Phase Eleven-Level Inverter              %
%         https://doi.org/10.1155/2020/8856040                            %
%_________________________________________________________________________%
clear
clc
SearchAgents_no=10; %  number of Black Widows Spiders 
Function_name='F28'; % <== write 'F2' or 'F3' and so on.  Name of the test function that can be from F1 to F23, for F24+ are ing. problems.
Max_iteration=100; % Maximum numbef of iterations
[lb,ub,dim,fobj]=Get_Function(Function_name);   % Load details of the selected benchmark function
[vMin,theBestVct,Convergence_curve]=BWOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
display(['The best solution obtained by BWOA is : ', num2str(theBestVct)]);
display(['The best optimal value of the objective funciton found by BWOA is : ', num2str(vMin)]);
%**********************************************************************************************************
        



