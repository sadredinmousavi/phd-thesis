function [Psai] = calculatePsai_3PM(point, MagPos, Psai1)
    
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
    Ar  = [
        c(2)*r_(2,1) c(3)*r_(3,1)
        c(2)*r_(2,2) c(3)*r_(3,2)
    ];
    br  = -[
        c(1)*r_(1,1)
        c(1)*r_(1,2)
    ];
    Psai1 = 30*pi/180;
    psai23 = acos(inv(Ar)*br*cos(Psai1))';
    psai123 = [Psai1 psai23];
    error = myFun(psai123');
    Psai = psai123;
    %
    
%     A = [];
%     b = [];
%     Aeq = [];
%     beq = [];
%     lb  = ones(1,size(MagPos,1))* 000 *pi/180;
%     ub  = ones(1,size(MagPos,1))* 090 *pi/180;
%     psai_0 = ones(1,size(MagPos,1))* 180 *pi/180;
%     psai_0 = lb;
%     [psai123,fval,exitflag,output] = fmincon(myFun2,psai_0',A,b,Aeq,beq,lb',ub');
%     error = myFun1(psai123);
%     Psai = psai123;

end