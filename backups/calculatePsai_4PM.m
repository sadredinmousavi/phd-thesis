function [Psai] = calculatePsai_4PM(point1, point2, MagPos)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(point1, MagPos);
    [r2, a2, b2, c2] = calculateParamsFromPoint(point2, MagPos);
    %
    %
    coeff = [r1'; r2'];
    myFun = @(psai)  norm( coeff*cos(psai) ) ;
    [V, D] = eig(coeff);
    Psai = acos( real( V(:,1) ));
    error = myFun(Psai);
    
    
    
    
    %
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb  = ones(1,size(MagPos,1))* 000 *pi/180;
    ub  = ones(1,size(MagPos,1))* 090 *pi/180;
    psai_0 = ones(1,size(MagPos,1))* 180 *pi/180;
    psai_0 = lb;
    [psai123,fval,exitflag,output] = fmincon(myFun,psai_0',A,b,Aeq,beq,lb',ub');
    error = myFun(psai123);
    Psai = psai123;

end