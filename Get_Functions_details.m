function [lb,ub,dim,fobj] = Get_Functions_details(F)
switch F       
    case 'F1'
        fobj = @F1;
        lb = [  1,    1,   -20,   -20,   -20,   -20,  -20,  -20,  -20 ];
        ub = [ 1e3,  1e3,   20,    20,    20,    20,   20,   20,   20 ]; 
        dim=9;  
end
end

function o = F1(x)
o=PI_RMRAC(x);
end


