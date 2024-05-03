function [r1, a1, b1, c1] = calculateParamsFromPoint(point, MagPos)

    x1 = point(1);
    y1 = point(2);
    for i=1:size(MagPos,1) % sum of all forces in point one
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r1_(i,:) = [x1-x_ y1-y_];
        r1__(i,1) = sqrt( r1_(i,1)^2 + r1_(i,2)^2 + (0-z_)^2 );
        r1i5(i,1) = 1/r1__(i)^5;
        r1i7(i,1) = 1/r1__(i)^7;
        r1(i,:) = r1_(i,:) * r1i5(i,1);
        a1(i,1) = ( -4*r1_(i,1)^2 + r1_(i,2)^2 ) * r1i7(i,1);
        b1(i,1) = ( -5*r1_(i,1) * r1_(i,2) ) * r1i7(i,1);
        c1(i,1) = ( +r1_(i,1)^2 - 4*r1_(i,2)^2 ) * r1i7(i,1);
    end
    
end

