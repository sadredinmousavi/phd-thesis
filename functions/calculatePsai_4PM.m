function [Psai] = calculatePsai_4PM(point1, point2, MagPos)
    
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
    coeff = [r1'; r2'];
    myFun = @(psai)  norm( coeff*cos(psai) ) ;
    [V, D] = eig(coeff);
    Psai = acos( real( V(:,1) ));
    error = myFun(Psai);
    
    
    
    
    %
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb  = ones(1,size(MagPos,1))* 000 *pi/180;
    ub  = ones(1,size(MagPos,1))* 090 *pi/180;
    psai_0 = ones(1,size(MagPos,1))* 180 *pi/180;
    psai_0 = lb;
    [psai123,fval,exitflag,output] = fmincon(myFun,psai_0',A,b,Aeq,beq,lb',ub');
    error = myFun(psai123);
    Psai = psai123;

end