function [rankM,error, hasAns, isStable, Psai, hessian, otherOutputs] = check_4PM(point1, point2, MagPos)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(point1, MagPos);
    [r2, a2, b2, c2] = calculateParamsFromPoint(point2, MagPos);
    %
    %
    coeff = [r1'; r2'];
    myFun = @(psai)  norm( coeff*cos(psai) ) ;
    [V, D] = eig(coeff);
    index = 4;
    for k=1:size(D,1)
        if D(k,k) == 0
            index = k;
        end
    end
    Psai = acos( V(:,index) );
    error = myFun(Psai);
    
    
%     A = [];
%     b = [];
%     Aeq = [];
%     beq = [];
%     lb  = ones(1,size(MagPos,1))* 000 *pi/180;
%     ub  = ones(1,size(MagPos,1))* 090 *pi/180;
%     psai_0 = ones(1,size(MagPos,1))* 45 *pi/180;
%     psai_0 = lb;
%     [psai123,fval,exitflag,output] = fmincon(myFun,psai_0',A,b,Aeq,beq,lb',ub');
%     error = myFun(psai123);
%     Psai = psai123;
    
    
    rankM = rank(coeff);
    if rankM < size(coeff,1)
        hasAns = 1;
    else
        hasAns = 0;
    end
%     error = myFun(psai123);
    
    [isStable1,hessian.point1] = isHessianStable(Psai, a1, b1, c1);
    [isStable2,hessian.point2] = isHessianStable(Psai, a2, b2, c2);
    if isStable1+isStable2 == 2
        isStable = 1;
    else
        isStable = 0;
    end
    
   otherOutputs = 0;

end