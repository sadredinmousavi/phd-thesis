function [r,isInContact,normalVector] = calculateDistanceToWallLinear(point, wall, threshold)
if nargin < 3
    threshold = 1000;
end
v1 = [ [wall(3:4,1) - wall(1:2,1)] ;0];  % vector from wall_1 to wall_2
v2 = [ [point(1:2,1) - wall(1:2,1)] ;0]; % vector from wall_1 to point
theta1 = atan2( norm(cross(v1,v2)), dot(v1,v2) ) *(180/pi);
%
v1 = [ [wall(1:2,1) - wall(3:4,1)] ;0];  % vector from wall_2 to wall_1
v2 = [ [point(1:2,1) - wall(3:4,1)] ;0]; % vector from wall_2 to point
theta2 = atan2( norm(cross(v1,v2)), dot(v1,v2) ) *(180/pi);
%
r = norm(v2)*sind(theta2);
if theta1 < 91 && theta2 < 91 && r < threshold
    isInContact = 1;
    projectedPointOnTheWall = wall(3:4,1) + v1(1:2)./norm(v1) * ( norm(v2)*sind(theta2) );
    normalVector = point(1:2,1) - projectedPointOnTheWall;
    normalVector = normalVector ./ norm(normalVector);
    % direction from wall to point
else
    isInContact = 0;
    normalVector = 0;
end

end
