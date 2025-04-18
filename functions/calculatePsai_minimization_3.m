function [rankM,error, hasAns, isStable, Psai, hessian, otherOutputs] = calculatePsai_minimization_3(points, MagPos, lambda, psai_0)
    
    
%     [r1, a1, b1, c1] = calculateParamsFromPoint(points{1}, MagPos);
%     [r2, a2, b2, c2] = calculateParamsFromPoint(points{2}, MagPos);
    
    [theta,rho] = cart2pol(points{1}(1), points{1}(2));
    point = [rho, 0];
    psai_num    = size(MagPos, 1);
    psai_r_num  = (psai_num/2 - 1) + 2;
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(point, MagPos);
    %
    coeff = r1';
    forceEqs = @(psai_r) coeff*cos([psai_r(1); psai_r(3:end); psai_r(2); psai_r(3:end)]);
    costFun = @(psai_r) norm( r1'*cos([psai_r(1); psai_r(3:end); psai_r(2); psai_r(3:end)]) );
    convertReducedToNormal = @(psai_r) [psai_r(1); psai_r(3:end); psai_r(2); psai_r(3:end)];
    %
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb  = [ -90; -90; -ones(psai_r_num-2,1)* 070 ] *pi/180 ;
    ub  = [ +90; +90; +ones(psai_r_num-2,1)* 070 ] *pi/180 ;
    Psai_0 = ones(psai_r_num,1)* 000 *pi/180 ;
%     Psai_0 = lb;
    %
    options = optimoptions('fmincon');
    options = optimoptions(options,'Display', 'Iter','Algorithm','interior-point');
    % [Psai1,fval,exitflag,output,lambda,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub,@(Psai) stabilityConstraints(Psai, points, MagPos), options);
    [Psai1_r,fval,exitflag,output,lambda,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub,@(Psai) stabilityConstraints(Psai, point, MagPos), options);
%     [Psai2,fval,exitflag,output,lambda,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub);

    Psai = convertReducedToNormal(Psai1_r);
    
    error = costFun(Psai1_r);
    error2 = forceEqs(Psai1_r);
    c_ = stabilityConstraints(Psai1_r, point, MagPos);

    rankM=4;
    hasAns=0;
    [isStable1,hessian1_] = isHessianStable_r(Psai1_r, a1, b1, c1);
    if isStable1 == 1
        isStable = 1;
    else
        isStable = 0;
    end
    otherOutputs = 0;
    
    function [c,ceq] = stabilityConstraints(Psai_r, point, MagPos)
        Psai_ = [Psai_r(1); Psai_r(3:end); Psai_r(2); Psai_r(3:end)];
        [r_1, a_1, b_1, c_1] = calculateParamsFromPoint(point, MagPos);
        hessian1 = [cos(Psai_)'*a_1 cos(Psai_)'*b_1;cos(Psai_)'*b_1 cos(Psai_)'*c_1];
        %
        [V,D] = eig(hessian1); % V(:,i)
        angle1(1,1) = atan(V(2,1)/V(1,1))*(180/pi);
        angle1(2,1) = atan(V(2,2)/V(1,2))*(180/pi);
        d1(1,1) = D(1,1);
        d1(2,1) = D(2,2);
        threshold = 1.3;
        if abs(d1(1)) > abs(d1(2))
            dist = ( abs( d1(1)/d1(2) ) - threshold ) * 10000;
        else
            dist = ( abs( d1(2)/d1(1) ) - threshold ) * 10000;
        end
        %
        %
        c(1:2) = [d1];
        c(3) = [dist];
        c(4) = costFun(Psai_r)*1000;
        ceq = [];
    end

    function [isStable,hessian1] = isHessianStable_r(Psai_r, a_1, b_1, c_1)
        Psai_ = [Psai_r(1); Psai_r(3:end); Psai_r(2); Psai_r(3:end)];
%         [r_1, a_1, b_1, c_1] = calculateParamsFromPoint(point, MagPos);
        hessian1 = [cos(Psai_)'*a_1 cos(Psai_)'*b_1;cos(Psai_)'*b_1 cos(Psai_)'*c_1];
        try 
            d = eig(hessian1);
        catch
            d = [1 1];
        end
        if d(1) < 0 && d(2) < 0 && ~all( round(Psai_,2) == round(pi/2,2) )
            isStable = 1;
        else
            isStable = 0;
        end
    end
    
end