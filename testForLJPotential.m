% sigma = 1.2 * 1e-5;
% epsilun = 0.01;
% LJ_potential = @(r)4*epsilun*( (sigma/r)^12 - (sigma/r)^6 );
% LJ_force     = @(r)4*epsilun*( -12*(sigma/r)^12/r + 6*(sigma/r)^6/r );
% % fplot(LJ_potential,[sigma/10,3*sigma]);
% % fplot(LJ_force,[0.000001,0.0005])
% fplot(LJ_force,[sigma/10,100*sigma])
% LJ_force(1.246*sigma)

wall = [-3;0;0;3];
point = [0;0];
[r,isInContact,normalVector] = calculateDistanceToWallLinear(point, wall)
plot(wall([1,3]), wall([2,4]))
hold on
plot(point(1), point(2), 'ro')