function [Psai] = calculatePsai_6PM(point1, point2, MagPos)
    
    x1 = point1(1);
    y1 = point1(2);
    x2 = point2(1);
    y2 = point2(2);
    
    if x1 == -1 && y1 == -1
        Psai = ones(1,size(MagPos,1))* 090 *pi/180;
        return;
    end

    for i=1:size(MagPos,1) % sum of all forces in point one
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r1_(i,:) = [x1-x_ y1-y_];
        r1__(i,1) = sqrt( r1_(i,1)^2 + r1_(i,2)^2 + (0-z_)^2 );
        c1(i,1) = 1/r1__(i)^5;
        r1(i,:) = r1_(i,:) * c1(i,1);
    end
    for i=1:size(MagPos,1) % sum of all forces in point two
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r2_(i,:) = [x2-x_ y2-y_];
        r2__(i,1) = sqrt( r2_(i,1)^2 + r2_(i,2)^2 + (0-z_)^2 );
        c2(i,1) = 1/r2__(i)^5;
        r2(i,:) = r1_(i,:) * c2(i,1);
    end
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
    ub  = ones(1,size(MagPos,1))* 090 *pi/180;
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

    
end