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
function [ val] = getBinary( )  % algorithm-2 in the paper page 6.
if rand() < 0.5
     val= 0;
else
     val=1;
end

