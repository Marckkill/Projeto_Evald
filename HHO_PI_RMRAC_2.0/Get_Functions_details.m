
function [lb,ub,dim,fobj] = Get_Functions_details(F)
switch F
    case 'F1'
        fobj = @F1;
        lb = [   1  500  -2 -2 -2 -2 -2 -2 -2 ];
        ub = [ 500 1000   2  2  2  2  2  2  2 ];
        dim= 9;             
end
end

% F1
function o = F1(x)

o = PI_RMRAC_for_optimization(x);

end
