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
function [ o ] =  getPheromone(  fit, min, max ) % Eq.12 in the paper
    for i=1:size(fit,2)
         o(i)= (max-fit(i))/(max-min);
    end
end

