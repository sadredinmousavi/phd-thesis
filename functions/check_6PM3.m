function [rankM,error, hasAns, isStable, Psai, hessian, otherOutputs] = check_6PM3(point1, point2, point3, MagPos)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(point1, MagPos);
    [r2, a2, b2, c2] = calculateParamsFromPoint(point2, MagPos);
    [r3, a3, b3, c3] = calculateParamsFromPoint(point3, MagPos);
    %
    %
    coeff = [r1'; r2'; r3'];
    myFun = @(psai)  norm( coeff*cos(psai) ) ;
    [V, D] = eig(coeff);
    Psai = acos( V(:,6) );
    error = myFun(Psai);
    
    
    rankM = rank(coeff);
    if rankM < size(coeff,1)
        hasAns = 1;
    else
        hasAns = 0;
    end
    
    [isStable1,hessian.point1] = isHessianStable(Psai, a1, b1, c1);
    [isStable2,hessian.point2] = isHessianStable(Psai, a2, b2, c2);
    [isStable3,hessian.point3] = isHessianStable(Psai, a3, b3, c3);
    
    if isStable1+isStable2+isStable3 == 3
        isStable = 1;
    else
        isStable = 0;
    end
    
    otherOutputs = 0;
end