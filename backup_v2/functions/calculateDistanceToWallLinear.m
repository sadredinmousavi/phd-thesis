function [r,isInContact,normalVector] = calculateDistanceToWallLinear(point, wall, threshold)
if nargin < 3
    threshold = 1000;
end
v1 = [ [wall(3:4,1) - wall(1:2,1)] ;0];  % vector from wall_1 to wall_2
v2 = [ [point(1:2,1) - wall(1:2,1)] ;0]; % vector from wall_1 to point
%
v3 = [ [wall(1:2,1) - wall(3:4,1)] ;0];  % vector from wall_2 to wall_1
v4 = [ [point(1:2,1) - wall(3:4,1)] ;0]; % vector from wall_2 to point
%
r = norm(cross(v1, v2)) / norm(v1);
if r < threshold
    theta1 = atan2( norm(cross(v1,v2)), dot(v1,v2) ) *(180/pi);
    theta2 = atan2( norm(cross(v3,v4)), dot(v3,v4) ) *(180/pi);
    %
    if theta1 <= 90 && theta2 <= 90
        isInContact = 1;
        projectedPointOnTheWall = wall(3:4,1) + v3(1:2)./norm(v3) * ( norm(v4)*cosd(theta2) );
        normalVector = point(1:2,1) - projectedPointOnTheWall;
        normalVector = normalVector ./ norm(normalVector);
        % direction from wall to point
    else
        isInContact = 0;
        normalVector = 0;
    end
else
    isInContact = 0;
    normalVector = 0;
end

end

