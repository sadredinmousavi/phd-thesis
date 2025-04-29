function [rankM,error, hasAns, isStable, Psai, hessian, otherOutputs] = calculatePsai_minimization_2(points, MagPos, lambda, psai_0)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(points{1}, MagPos);
    coeff = r1';
    for i=2:length(points)
        [r1, a1, b1, c1] = calculateParamsFromPoint(points{i}, MagPos);
        coeff = [coeff; r1'];
    end
    %
    forceEqs = @(psai) coeff*cos(psai);
    costFun = @(psai) norm( coeff*cos(psai) )  + max(mean(abs(psai)) - 35*(pi/180), 0) * 1e5;
%     costFun = @(psai) norm( coeff*cos(psai) );
%     costFun = @(psai) norm( r1'*cos(psai) ) + norm( r2'*cos(psai) ) + norm( r3'*cos(psai) );


    %
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb  = -ones(size(MagPos,1),1)* 090 *pi/180 ;
    ub  = +ones(size(MagPos,1),1)* 090 *pi/180 ;
    if nargin > 3
        Psai_0 = psai_0;
%         Psai_0 = [-85.00 -33.00 -33.00 -85.00 -33.00 -33.00]'*(pi/180);
%         lb  = 0.90*Psai_0;
%         ub  = 1.10*Psai_0;
    else
        Psai_0 = 0.1*lb;
    end
    %
    options = optimoptions('fmincon');
    options = optimoptions(options,'Display', 'off','Algorithm','sqp', 'MaxFunctionEvaluations', 5e3); % 
    if nargin > 2
        [Psai1,fval,exitflag,output,lambda1,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub,@(Psai) stabilityConstraints(Psai, points, MagPos, lambda), options);
    else
        [Psai1,fval,exitflag,output,lambda1,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub,@(Psai) stabilityConstraints(Psai, points, MagPos), options);
%         [Psai2,fval,exitflag,output,lambda,grad,hessian] = fmincon(costFun,Psai_0,A,b,Aeq,beq,lb,ub);
    end
    

    Psai = Psai1;
    
    error = costFun(Psai);
    error2 = forceEqs(Psai);
    if nargin > 2
        c_ = stabilityConstraints(Psai, points, MagPos, lambda);
    else
        c_ = stabilityConstraints(Psai, points, MagPos);
    end

    rankM=4;
    hasAns=0;
    stableNum = 0;
    for i=1:length(points)
        [r1, a1, b1, c1] = calculateParamsFromPoint(points{i}, MagPos);
        [isStable1,hessian1_] = isHessianStable(Psai, a1, b1, c1);
        if isStable1
            stableNum = stableNum + 1;
        end
    end    
    if stableNum == length(points)
        isStable = 1;
    else
        isStable = 0;
    end
    otherOutputs = 0;
    
    function [c,ceq] = stabilityConstraints(Psai, points, MagPos, extraParam)
        c = zeros(2*length(points), 1);
        ceq = [];
        for k=1:length(points)
            [r_1, a_1, b_1, c_1] = calculateParamsFromPoint(points{k}, MagPos);
            hessian1 = [cos(Psai)'*a_1 cos(Psai)'*b_1;cos(Psai)'*b_1 cos(Psai)'*c_1];
            %
            [V,D] = eig(hessian1); % V(:,i)
            angle1(1,1) = atan(V(2,1)/V(1,1))*(180/pi);
            angle1(2,1) = atan(V(2,2)/V(1,2))*(180/pi);
            d1(1,1) = D(1,1);
            d1(2,1) = D(2,2);
            c(2*k-1:2*k, 1) = d1;
        end
        if nargin > 3
            ref = angle1(1,1);
            const = abs( ref - extraParam{1} );
            if const < 0.1 * ref
                const = - const;
            end
            c = [c; angle1(1,1)-extraParam{1}];
        end
%         k=k+1;
%         [r_1, a_1, b_1, c_1] = calculateParamsFromPoint([0,0], MagPos);
%         hessian1 = [cos(Psai)'*a_1 cos(Psai)'*b_1;cos(Psai)'*b_1 cos(Psai)'*c_1];
%         %
%         [V,D] = eig(hessian1); % V(:,i)
%         angle1(1,1) = atan(V(2,1)/V(1,1))*(180/pi);
%         angle1(2,1) = atan(V(2,2)/V(1,2))*(180/pi);
%         d1(1,1) = -D(1,1);
%         d1(2,1) = -D(2,2);
%         c(2*k-1:2*k, 1) = d1;
    end
    
end