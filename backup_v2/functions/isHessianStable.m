function [isStable,hessian] = isHessianStable(Psai, a, b, c)
    hessian = [cos(Psai)'*a cos(Psai)'*b;cos(Psai)'*b cos(Psai)'*c];
    try 
        d = eig(hessian);
    catch
        d = [1 1];
    end
    if d(1) < 0 && d(2) < 0 && ~all( round(Psai,2) == round(pi/2,2) )
        isStable = 1;
    else
        isStable = 0;
    end
%     if d(1) > 0 && d(2) > 0 && ~all( round(Psai,2) == round(pi/2,2) )
%         isStable = 1;
%     end
end