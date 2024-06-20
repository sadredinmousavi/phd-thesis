function [rankM,error, hasAns, isStable, Psai, hessian, otherOutputs] = calculatePsai_minimization(points, MagPos, lambda, psai_0)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(points{1}, MagPos);
    [r2, a2, b2, c2] = calculateParamsFromPoint(points{2}, MagPos);
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
    lb  = -ones(size(MagPos,1),1)* 090 *pi/180 ;
    ub  = +ones(size(MagPos,1),1)* 070 *pi/180 ;
%     Psai_0 = ones(size(MagPos,1),1)* 180 *pi/180 ;
    Psai_0 = lb*0.001;
    %
    options = optimoptions('fmincon');
    options = optimoptions(options,'Display','off','Algorithm','interior-point'); % Iter
    [Psai1,fval,exitflag,output,lambda,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub,@(Psai) stabilityConstraints(Psai, points, MagPos), options);
%     [Psai2,fval,exitflag,output,lambda,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub);

    Psai = Psai1;
    
    error = costFun(Psai);
    error2 = forceEqs(Psai);
    c_ = stabilityConstraints(Psai, points, MagPos);

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
    
    function [c,ceq] = stabilityConstraints(Psai, points, MagPos)
        [r_1, a_1, b_1, c_1] = calculateParamsFromPoint(points{1}, MagPos);
        hessian1 = [cos(Psai)'*a_1 cos(Psai)'*b_1;cos(Psai)'*b_1 cos(Psai)'*c_1];
        [r_2, a_2, b_2, c_2] = calculateParamsFromPoint(points{2}, MagPos);
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
        threshold = 1.3;
        if abs(d1(1)) > abs(d1(2))
            dist1 = ( abs( d1(1)/d1(2) ) - threshold ) * 100;
        else
            dist1 = ( abs( d1(2)/d1(1) ) - threshold ) * 100;
        end
        if abs(d2(1)) > abs(d2(2))
            dist2 = ( abs( d2(1)/d2(2) ) - threshold ) * 100;
        else
            dist2 = ( abs( d2(2)/d2(1) ) - threshold ) * 100;
        end
        %
        c(1:4) = [d1; d2];
        c(5)   = costFun(Psai);
        if points{1}(1) == points{2}(1) && points{1}(2) == points{2}(2)
            c(6:7) = [dist1; dist2];
        end
        ceq = [];
    end
    
end