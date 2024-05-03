function [rankM,error, hasAns, isStable, Psai, hessian, otherOutputs] = calculatePsai_6PM(points, MagPos, lambda, psai_0)
    
    [r1, a1, b1, c1] = calculateParamsFromPoint(points{1}, MagPos);
    
    [r2, a2, b2, c2] = calculateParamsFromPoint(points{2}, MagPos);
    %
    %
    coeff = [r1'; r2'];
    myFun = @(psai)  norm( coeff*cos(psai) ) ;
    myFun2 = @(psai) norm( r1'*cos(psai) ) + norm( r2'*cos(psai) );
      
    index1 = 2; % indexes of two controlling actuators which goes other side of equality
    index2 = 6;
    ind1 = setdiff(1:size(coeff,2),[index1, index2]);
    ind2 = [index1, index2];
    Ar = coeff(:,ind1);
    Br = -coeff(:,ind2);
    
    hasAns = 0;
    ansCount = 0;
    
    for psai1=(-90:1:90)*(pi/180)
        for psai2=(-90:1:90)*(pi/180)
            psai_temp = acos( inv(Ar)*Br*cos([psai1;psai2]) );
            psai_selc = [psai1;psai2];
            if isreal(psai_temp) && rank(Ar) == 4
                hasAns = 1;
                for cnt=1:size(coeff,2)
                    if any(ind1(:) == cnt) %is inside Ar
                        Psai(cnt,1) = psai_temp(find(ind1(:) == cnt));
                    else % is inside psai1 and psai 2
                        Psai(cnt,1) = psai_selc(find(ind2(:) == cnt));
                    end
                end
                [isStable1,hessian.point1] = isHessianStable(Psai, a1, b1, c1);
                [isStable2,hessian.point2] = isHessianStable(Psai, a2, b2, c2);
                if isStable1+isStable2 == 2
                    isStable = 1;
                else
                    isStable = 0;
                    continue
                end
                
                break
                ansCount = ansCount + 1;
                PsaiSerie(:,ansCount) = Psai;
                [V,D] = eig(hessian.point1); % V(:,i)
                angle1_point1(ansCount) = atan(V(2,1)/V(1,1))*(180/pi);
                angle2_point1(ansCount) = atan(V(2,2)/V(1,2))*(180/pi);
                d1_point1(ansCount) = D(1,1);
                d2_point1(ansCount) = D(2,2);
                [V,D] = eig(hessian.point2); % V(:,i)
                angle1_point2(ansCount) = atan(V(2,1)/V(1,1))*(180/pi);
                angle2_point2(ansCount) = atan(V(2,2)/V(1,2))*(180/pi);
                d1_point2(ansCount) = D(1,1);
                d2_point2(ansCount) = D(2,2);
            else
                hasAns = 0;
            end
        end
        if hasAns && isStable
%         if hasAns
            break
        end
    end
    
    otherOutputs = 0;
%     otherOutputs.PsaiSerie = PsaiSerie;
%     otherOutputs.angle1_point1 = angle1_point1;
%     otherOutputs.angle2_point1 = angle2_point1;
%     otherOutputs.d1_point1 = d1_point1;
%     otherOutputs.d2_point1 = d2_point1;
%     otherOutputs.angle1_point2 = angle1_point2;
%     otherOutputs.angle2_point2 = angle2_point2;
%     otherOutputs.d1_point2 = d1_point2;
%     otherOutputs.d2_point2 = d2_point2;
    
    rankM = rank(Ar);
    error = myFun(Psai);
    
    [isStable1,hessian.point1] = isHessianStable(Psai, a1, b1, c1);
    [isStable2,hessian.point2] = isHessianStable(Psai, a2, b2, c2);
    if isStable1+isStable2 == 2
        isStable = 1;
    else
        isStable = 0;
    end
    
end