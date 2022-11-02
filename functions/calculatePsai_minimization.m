function [rankM,error, hasAns, isStable, Psai, hessian, otherOutputs] = calculatePsai_minimization(point1, point2, MagPos)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(point1, MagPos);
    [r2, a2, b2, c2] = calculateParamsFromPoint(point2, MagPos);
    %
    coeff = [r1'; r2'];
    forceEqs = @(psai) coeff*cos(psai);
%     costFun = @(psai) norm( coeff*cos(psai) );
    costFun = @(psai) norm( r1'*cos(psai) ) + norm( r2'*cos(psai) );
    %
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb  = ones(size(MagPos,1),1)* 000 *pi/180 ;
    ub  = ones(size(MagPos,1),1)* 180 *pi/180 ;
%     Psai_0 = ones(size(MagPos,1),1)* 180 *pi/180 ;
    Psai_0 = lb;
    %
    options = optimoptions('fmincon');
    options = optimoptions(options,'Display', 'off');
    [Psai1,fval,exitflag,output,lambda,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub,@(Psai) stabilityConstraints(Psai, point1, point2, MagPos), options);
    [Psai2,fval,exitflag,output,lambda,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub);

    Psai = Psai1;
    
    error = costFun(Psai);
    error2 = forceEqs(Psai);
    c_ = stabilityConstraints(Psai, point1, point2, MagPos);

    rankM=4;
    hasAns=0;
    [isStable1,hessian1_] = isHessianStable(Psai, a1, b1, c1);
    [isStable2,hessian2_] = isHessianStable(Psai, a2, b2, c2);
    if isStable1+isStable2 == 2
        isStable = 1;
    else
        isStable = 0;
    end
    otherOutputs = 0;
    
    function [c,ceq] = stabilityConstraints(Psai, point1, point2, MagPos)
        [r_1, a_1, b_1, c_1] = calculateParamsFromPoint(point1, MagPos);
        hessian1 = [cos(Psai)'*a_1 cos(Psai)'*b_1;cos(Psai)'*b_1 cos(Psai)'*c_1];
        [r_2, a_2, b_2, c_2] = calculateParamsFromPoint(point2, MagPos);
        hessian2 = [cos(Psai)'*a_2 cos(Psai)'*b_2;cos(Psai)'*b_2 cos(Psai)'*c_2];
        %
        [V,D] = eig(hessian1); % V(:,i)
        angle1(1,1) = atan(V(2,1)/V(1,1))*(180/pi);
        angle1(2,1) = atan(V(2,2)/V(1,2))*(180/pi);
        d1(1,1) = D(1,1);
        d1(2,1) = D(2,2);
        %
        [V,D] = eig(hessian2); % V(:,i)
        angle2(1,1) = atan(V(2,1)/V(1,1))*(180/pi);
        angle2(2,1) = atan(V(2,2)/V(1,2))*(180/pi);
        d2(1,1) = D(1,1);
        d2(2,1) = D(2,2);
        %
        c(1:4) = [d1; d2];
        ceq = [];
    end
    
end