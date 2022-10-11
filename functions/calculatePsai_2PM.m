function [Psai] = calculatePsai_2PM(point, MagPos)
    
    x1 = point(1);
    y1 = point(2);
    
    if x1 == -1 && y1 == -1
        Psai = ones(1,size(MagPos,1))* 090 *pi/180;
        return;
    end

    for i=1:size(MagPos,1)
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r_(i,:) = [x1-x_ y1-y_];
        r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
        c(i,1) = 1/r__(i)^5;
        r(i,:) = r_(i,:) * c(i,1);
    end
    %
    %
    myFun = @(psai) norm( r'*cos(psai) );
    coeff = r';
    [V, D] = eig(coeff);
    Psai = acos(V(:,2));
    error = myFun(Psai);
    %
%     A = [];
%     b = [];
%     Aeq = [];
%     beq = [];
%     lb  = ones(1,size(MagPos,1))* 000 *pi/180;
%     ub  = ones(1,size(MagPos,1))* 090 *pi/180;
%     psai_0 = lb;
%     [psai123,fval,exitflag,output] = fmincon(myFun,psai_0',A,b,Aeq,beq,lb',ub');
%     error = myFun(psai123);
%     Psai = psai123;
    
end