function [lb,ub,dim,fobj] = Get_Functions_details(F)
switch F
    case 'F1'
        fobj = @F1;
        % Kp, Ki, Gamma, Kappa, Thetac(4), Thetas(4)
        lb= [ 1  1  1    1  -1 -1 ]; 
        ub= [ 20   20   100 2000  40  40 ];
        dim=6;    
        
        case 'F2'
        fobj = @F2;
        % Kp, Ki, Kd
        lb= [  1  1  0.0001     ];
        ub= [ 20 20 1  ];
        dim=3; 
end
end

% F1
function o = F1(x)

o = LLCL_for_opt(x);

end

function o = F2(x)

o = plant(x);

end
