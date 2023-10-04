%  Please refer to the main paper:
% Ali Asghar Heidari, Seyedali Mirjalili, Hossam Faris, Ibrahim Aljarah, Majdi Mafarja, Huiling Chen
% Harris hawks optimization: Algorithm and applications
% Future Generation Computer Systems, DOI: https://doi.org/10.1016/j.future.2019.02.028
% _______
function [lb,ub,dim,fobj] = Get_Functions_details(F)
switch F
    case 'F1'
        fobj = @F1;
        lb= [ 1 1 0.0001 ];
        ub= [ 20 20 20 ];
        dim=3;             
end
end

% F1
function o = F1(x)

o = plant(x);

end
