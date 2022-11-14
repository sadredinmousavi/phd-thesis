function [rankM,error, hasAns, isStable, Psai, hessian, otherOutputs] = calculatePsai_6PM(points, MagPos)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(points{1}, MagPos);
    [r2, a2, b2, c2] = calculateParamsFromPoint(points{2}, MagPos);
    %
    %
    myFun2 = @(psai) norm( r1'*cos(psai) ) + norm( r2'*cos(psai) );
    psai_0 = [30 150 181]*pi/180;
%     options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
%     [psai123,fval,exitflag,output] = fsolve(myFun1,psai_0,options);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb  = ones(1,size(MagPos,1))* 000 *pi/180;
    ub  = ones(1,size(MagPos,1))* 180 *pi/180;
    psai_0 = ones(1,size(MagPos,1))* 180 *pi/180;
    psai_0 = lb;
%     [psai123,fval,exitflag,output] = fminsearch(@(psai) abs(myFun2(psai)),psai_0);
    [psai123,fval,exitflag,output] = fmincon(myFun2,psai_0',A,b,Aeq,beq,lb',ub');
%     error = myFun1(psai123);
    error = myFun2(psai123);
    Psai = psai123;
    
    %
    %
%     Ar  = [
%         c3*r_(3,1) c2*r_(2,1)
%         c3*r_(3,2) c2*r_(2,2)
%     ];
%     br  = [
%         c1*r_(1,1)
%         c1*r_(1,2)
%     ];
%     psai23 = acos(inv(Ar)*br*cos(Psai1))';
%     psai123 = [Psai1 psai23];
%     error = myFun1(psai123);
% %     Psai = psai123;
%     %
    
    rankM=4;
    hasAns=0;
    [isStable1,hessian.point1] = isHessianStable(Psai, a1, b1, c1);
    [isStable2,hessian.point2] = isHessianStable(Psai, a2, b2, c2);
    if isStable1+isStable2 == 2
        isStable = 1;
    else
        isStable = 0;
    end
    
    
    otherOutputs = 0;
    
end