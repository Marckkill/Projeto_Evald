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

function [vMin,theBestVct,Convergence_curve]=BWOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)
Positions=initialization(SearchAgents_no,dim,ub,lb);
 for i=1:size(Positions,1)
      Fitness(i)=fobj(Positions(i,:)); % get fitness     
 end
[vMin minIdx]= min(Fitness);  % the min fitness value vMin and the position minIdx
theBestVct= Positions(minIdx,:);  % the best vector 
[vMax maxIdx]= max(Fitness);  % the min fitness value vMin and the position minIdx
Convergence_curve=zeros(1,Max_iter);
Convergence_curve(1)= vMin;
pheromone= getPheromone( Fitness, vMin, vMax);
t=0;% Loop counter
% Main loop
for t=1:Max_iter  
    beta= -1 + 2* rand();  % -1 < beta2 < 1     section 3.2.1 in the paper, page 4
    m= 0.4 + 0.5 *rand();  % 0.4 < m < 0.9 
   for r=1:SearchAgents_no       
      P= rand();
       r1= round(1+ (SearchAgents_no-1)* rand());
     if P >= 0.3  % spiral search   Eq. 11 in the paper, page 4
         v(r,:)=   theBestVct -  cos(2*pi*beta)*Positions(r,:);
     else         % direct search Eq. 11
          v(r,:)=   theBestVct - m *Positions(r1,:);
     end
    if pheromone(r) <= 0.3
         band=1; 
         while band 
           r1= round(1+ (SearchAgents_no-1)* rand());
           r2= round(1+ (SearchAgents_no-1)* rand());
           if r1 ~= r2 
               band=0;
           end
         end             % pheromone function.  Eq. 13 page 5 , getBinary is the  algorithm-2 in the paper, page 6
              v(r,:)=   theBestVct + (Positions(r1,:)-((-1)^getBinary)*Positions(r2,:))/2;
    end    
   %********************************************************************************
     % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=v(r,:)>ub;
        Flag4lb=v(r,:)<lb;
        v(r,:)=(v(r,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    % Evaluate new solutions
     Fnew= fobj(v(r,:));
     % Update if the solution improves
     if Fnew <= Fitness(r);
        Positions(r,:)= v(r,:);
        Fitness(r)= Fnew;
     end
     if Fnew <= vMin
         theBestVct= v(r,:);
         vMin= Fnew;
     end    
   end
    % update max and pheromons
   [vMax maxIdx]= max(Fitness);
   pheromone= getPheromone( Fitness, vMin, vMax);
   Convergence_curve(t+1)= vMin; 
 end
%**********************************[End  BWOA function]


