function [Psai] = calculatePsai_3PM(x, y, MagPos, Psai1)

    for i=1:size(MagPos,1)
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r_(i,:) = [x-x_ y-y_];
        r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
    end
    %
    %
    c1 = 1/r__(1)^5;
    c2 = 1/r__(2)^5;
    c3 = 1/r__(3)^5;
    myFun1 = @(psai)[c1*r_(1,1)*cos(psai(1))+c2*r_(2,1)*cos(psai(2))+c3*r_(3,1)*cos(psai(3)) c1*r_(1,2)*cos(psai(1))+c2*r_(2,2)*cos(psai(2))+c3*r_(3,2)*cos(psai(3))];
    myFun2 = @(psai) sqrt( (c1*r_(1,1)*cos(psai(1))+c2*r_(2,1)*cos(psai(2))+c3*r_(3,1)*cos(psai(3)))^2 +  (c1*r_(1,2)*cos(psai(1))+c2*r_(2,2)*cos(psai(2))+c3*r_(3,2)*cos(psai(3)))^2 );
    psai_0 = [30 01 181]*pi/180;
%     options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
%     [psai123,fval,exitflag,output] = fsolve(myFun1,psai_0,options);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb  = [0,0,0]*pi/180;
    ub  = [360,360,360]*pi/180;
%     [psai123,fval,exitflag,output] = fminsearch(@(psai) abs(myFun2(psai)),psai_0);
    [psai123,fval,exitflag,output] = fmincon(myFun2,psai_0,A,b,Aeq,beq,lb,ub);
    error = myFun1(psai123);
    Psai = psai123;
    %
    %
    Ar  = [
        c3*r_(3,1) c2*r_(2,1)
        c3*r_(3,2) c2*r_(2,2)
    ];
    br  = [
        c1*r_(1,1)
        c1*r_(1,2)
    ];
    psai23 = acos(inv(Ar)*br*cos(Psai1))';
    psai123 = [Psai1 psai23];
    error = myFun1(psai123);
%     Psai = psai123;
    %

    
end