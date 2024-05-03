function [r,theta] = calculateDistanceToWall(point, wall)
v1 = wall(1:2,1) - wall(3:4,1);
v2 = wall(1:2,1) - point(1:2,1);
theta1 = atan2( norm(cross(v1,v2)), dot(v1,v2) );
%
v1 = wall(1:2,1) - wall(3:4,1);
v2 = wall(3:4,1) - point(1:2,1);
theta2 = atan2( norm(cross(v1,v2)), dot(v1,v2) );
%
r=0;
theta=0;
end

