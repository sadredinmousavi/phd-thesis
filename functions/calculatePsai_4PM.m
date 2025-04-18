function [rankM, error, hasAns, isStable, Psai, hessian, otherOutputs] = calculatePsai_4PM(points, MagPos, lambda, psai_0)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(points{1}, MagPos);
    coeff = r1';
    for i=2:length(points)
        [r1, a1, b1, c1] = calculateParamsFromPoint(points{i}, MagPos);
        coeff = [coeff; r1'];
    end
    %
    %
    myFun = @(psai)  norm( coeff*cos(psai) ) ;
    [V, D] = eig(coeff);
    indice = 1;
    for i=1:size(coeff, 1)
        if D(i,i) == 0
            indice = i;
        end
    end
    Psai = acos( real( V(:,indice) ));
    error = myFun(Psai);
    
    
    
    
    %
    
%     A = [];
%     b = [];
%     Aeq = [];
%     beq = [];
%     lb  = ones(1,size(MagPos,1))* 000 *pi/180;
%     ub  = ones(1,size(MagPos,1))* 090 *pi/180;
%     lb  = - ub;
%     psai_0 = ones(1,size(MagPos,1))* 030 *pi/180;
% %     psai_0 = lb;
%     [psai123,fval,exitflag,output] = fmincon(myFun,psai_0',A,b,Aeq,beq,lb',ub');
%     error = myFun(psai123);
%     Psai = psai123;
    
    %
    hasAns = 1;
    isStable = 1;
    hessian = 1;
    rankM = 1;
    otherOutputs = 1;

end