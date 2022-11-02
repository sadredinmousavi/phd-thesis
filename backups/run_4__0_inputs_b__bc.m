%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 
 % field lines of an axially magnetized cylinder. 

 
 % run symboic.m and f -> force_field_symbolic.m
 
 % todo MagPosR -> input to symbolic.m ; automatically generate force_field_symbolic.m
clc, close all, clearvars -except cursor_info
%%
mu_0 = 4*pi*1e-7;
M = 1.2706/mu_0;              % Magnetization   [A/m]
% magnets
L = 0.004;
D = 0.002;
m = M * (pi*D^2/4*L);
args.L = L;
args.D = D;
% agents
L = 0.0004;
D = 0.0002;
m_agents = M * (pi*D^2/4*L);
m_agents = m_agents * 1e3; %%%%%%%%% note
%
%

%[x y z m_norm mu_0 m_agents]
% values = [
% %     0.06 0.06 0.1 m mu_0 m_agents
% %     -0.06 0.06 0.1 m mu_0 m_agents
% %     -0.06 -0.06 0.1 m mu_0 m_agents
% %     0.06 -0.06 0.1 m mu_0 m_agents
% %     0.06*cosd(30) 0.06*sind(30) 0.1 m mu_0 m_agents
%     -0.06*cosd(30) 0.06*sind(30) 0.1 m mu_0 m_agents
%     0 -0.09 0.1 m mu_0 m_agents
%     
% ];
a = 0.25;
z = 0.1;
values = [
    a*cosd(30) a*sind(30) z m mu_0 m_agents
    -a*cosd(30) a*sind(30) z m mu_0 m_agents
    0 -a z m mu_0 m_agents
];
% values = [
%     0 +a z m mu_0 m_agents
%     0 -a z m mu_0 m_agents
% ];


MagPos = [values(:,1:3) zeros(size(values,1),3)];
for i=1:size(MagPos, 1)
    MagPos(i,4:6) = round(MagPos(i,4:6) ./ norm(MagPos(i,4:6)), 5);
end


Psai1 = 40*pi/180;
Point = [-0.00 0.01];
%%%
% Psai_0 = calcPsai2(Point(1), Point(2), MagPos);
%%%
Psai_0 = calcPsai3(Point(1), Point(2), MagPos, Psai1);
%%%


% Psai_0 = zeros(size(values,1),1);
% Psai_0 = [
% %     30
%     150
%     -90
% ]*pi/180;
% Psai_0 = Psai_0 + [
% %     0
%     0
%     -75
% ]*pi/180;



%
%

spaceRegion = -0.3:0.01:0.3;
x = spaceRegion; z = x; y = x;
Npoints = length(x);

%
%

tspan = 0:1/5:30;%tspan = 0:1/24:30;
% x = [ -0.1 0.05 0.1];
% y = [ 0.1 0.05 0.15];
% x_ = [ 0.06 0.05 0.04 0.061 0.051 0.041];
% y_ = [ 0.06 0.05 0.04 0.061 0.051 0.041];
x__ = min(x)+0.1:0.1:max(x)-0.1;
y__ = min(y)+0.1:0.1:max(y)-0.1;
[x_,y_] = meshgrid(x__,y__);
x_ = reshape(x_, 1, []);
y_ = reshape(y_, 1, []);

%
%
counterName = 1;
cd('data')
fileName1 = '';
while strcmp(fileName1, '')
    if isfile([num2str(counterName) '_datas.mat'])
        counterName = counterName + 1;
    else
        fileName1 = [num2str(counterName) '_input.mat'];
        fileName2 = [num2str(counterName) '_movie.gif'];
        fileName3 = [num2str(counterName) '_datas.mat'];
        fileName4 = [num2str(counterName) '_cursor.mat'];
        fileName5 = [num2str(counterName) '_all_vars.mat'];
    end
end
cd('..')



dlgTitle    = 'PreCheck';
dlgQuestion = 'Do you want to pre check?';
choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
if strcmp(choice, 'Yes') == 1
    DoHandPreCheck = 1;
    symbolic_fun_b(values, fileName3);
else
    DoHandPreCheck = 0;
end


[Fx, Fy, Fz] = CylFfield3(Psai_0,Point);
[Bx, By, Bz] = CylBfield3(Psai_0,Point);


[F, Frho, Faxial] = calculateMagneticForce(x,y,Psai_0);
[B, Brho, Baxial] = calculateMagneticField(x,y,Psai_0);

[eq_x,eq_y] = findEqPoints(x,y,Psai_0);




%%



save('input_data')
cd('data')
save(fileName5)
cd('..')



% printFig('input_data', Psai_0, [num2str(counterName) 'a']);


% run_4__4
% run_4__3
if ~DoHandPreCheck
%     run_4__3
end




%%





function [Bx, By, Bz] = CylBfield3(Psai,Point)
    B_ = magnetic_field_symbolic(Point(1), Point(2), Psai);
    Bx = B_(1);
    By = B_(2);
    Bz = B_(3);
end
%
function [Fx, Fy, Fz] = CylFfield3(Psai,Point)
    F = force_field_symbolic(Point(1), Point(2), Psai);
    Fx = F(1);
    Fy = F(2);
    Fz = F(3);
end
%
function [B, Brho, Baxial] = calculateMagneticField(x,y,Psai)
    Npoints = length(x);
    Brho  = zeros(Npoints,Npoints); Baxial  = Brho; B  = Brho;
    %
    for i = 1:Npoints
        for j = 1:Npoints
            Point = [x(j) y(i) 0];
            [Bx1, By1, ~] = CylBfield3(Psai,Point');
            %
            Brho(i,j)   = Bx1 + Brho(i,j);
            Baxial(i,j) = By1 + Baxial(i,j);
        end
    end
    B = sqrt(Brho.^2+Baxial.^2);
end
%
function [F, Frho, Faxial] = calculateMagneticForce(x,y,Psai)
    Npoints = length(x);
    Frho = zeros(Npoints,Npoints); Faxial = Frho; F1 = Frho;
    %
    for i = 1:Npoints
        for j = 1:Npoints
            Point = [x(j) y(i) 0];
            [Fx1, Fy1, ~] = CylFfield3(Psai,Point');
            %
            Frho(i,j)   = Fx1 + Frho(i,j);
            Faxial(i,j) = Fy1 + Faxial(i,j);
        end
    end
    F = sqrt(Frho.^2+Faxial.^2);
end
%
function [eq_x,eq_y] = findEqPoints(x, y, Psai, maxTime) % Frho ro begire va ba delF hesab kone na inke dynamics bere psiCont ro brgire
    
    if nargin < 4
        maxTime = 150;
    end
    
    x__ = min(x)+0.0:0.1:max(x)-0.0;
    y__ = min(y)+0.0:0.1:max(y)-0.0;
    [x_,y_] = meshgrid(x__,y__);
    x_ = reshape(x_, 1, []);
    y_ = reshape(y_, 1, []);
%     x_ = 0.1;
%     y_ = 0.1;
    %
    ans0 = [x_ y_ zeros(1,length(x_)) zeros(1,length(x_))];
    [t,ans1] = ode45(@(t,y) systemDynamicsSimplified(t,y,Psai), [0 maxTime], ans0);
    n = size(ans1,2)/4;
    %
%     ans0 = [x_ y_];
%     ans1 = ans0;
%     n = length(x_);
%     [Frho_, Faxial_, ~] = CylFfield3(Psai,[ans0(1),ans1(1) 0]');
%     gama = ones(1,n)/Frho_/10000;
%     for i=1:100
%         for j = 1:n
%             [Frho(1,j), Faxial(1,j), ~] = CylFfield3(Psai,[ans1(j),ans1(n+j) 0]');
%         end
%         ans1(i+1,:) = ans1(i,:) + [gama.*Frho gama.*Faxial];
%     end
%     
    %
    [~, Frho, Faxial] = calculateMagneticForce(x,y,Psai);
    plot(ans1(end,1:n),ans1(end,n+1:2*n),'b.','MarkerSize', 16);
    hold on
    for j=1:n
        plot(ans1(:,j),ans1(:,n+j),'r-', 'LineWidth',1 );
    end
    plot_field = streamslice(x,y,Frho,Faxial,'method','cubic');
    set(plot_field,'Color','black','LineWidth',1.2);
    %
    counter = 0;
    for j=1:n
        if ans1(end,j) == ans1(end-1,j)
            plot(ans1(end,j),ans1(end,n+j),'g.','MarkerSize', 16);
            counter = counter + 1;
            eq_x(counter) = ans1(end,j); 
            eq_y(counter) = ans1(end,n+j);
        end
    end
    [C, ia, ic] = unique(round(eq_x,5));
    eq_x = eq_x(ia);
    eq_y = eq_y(ia);
    %
    a=1;
end
%
function [MagPos] = sphericalToCartesian(MagPosR)

    for i_=1:size(MagPosR,1)
        r = MagPosR(i_,1);
        theta = MagPosR(i_,2);
        phi = MagPosR(i_,3);
        psai = MagPosR(i_,4);

        x_ = r*sin(theta)*cos(phi);
        y_ = r*sin(theta)*sin(phi);
        z_ = r*cos(theta);
        u_ = [x_ y_ z_]/r;

        x__ = -cos(theta)*cos(phi);
        y__ = -cos(theta)*sin(phi);
        z__ = sin(theta);
        u__ = [x__ y__ z__];

        R_mat = [
            cos(psai)+u_(1)^2*(1-cos(psai))  u_(1)*u_(2)*(1-cos(psai))-u_(3)*sin(psai) u_(1)*u_(3)*(1-cos(psai))+u_(2)*sin(psai);
            u_(1)*u_(2)*(1-cos(psai))+u_(3)*sin(psai)  cos(psai)+u_(2)^2*(1-cos(psai)) u_(2)*u_(3)*(1-cos(psai))-u_(1)*sin(psai);
            u_(1)*u_(3)*(1-cos(psai))-u_(2)*sin(psai) u_(2)*u_(3)*(1-cos(psai))+u_(1)*sin(psai) cos(psai)+u_(3)^2*(1-cos(psai));
        ];

        MagPos(i_,:) = [x_ y_ z_ (R_mat*u__')'];
    end

end
%
function [Psai] = calcPsai2(x, y, MagPos)

    for i=1:size(MagPos,1)
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r_(i,:) = [x-x_ y-y_];
        r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
    end
    c = -r__(1)^5/r__(2)^5;
    myFun = @(psai2)asin(c*sin(psai2))-acos(c*cos(psai2));
    psai2_0 = 40*pi/180;
    [psai2,fval,exitflag,output] = fsolve(myFun,psai2_0);
    [psai2,fval,exitflag,output] = fminsearch(@(psai2) abs(myFun(psai2)),psai2_0);
    psai1 = asin(c*sin(psai2));
    psai1 = acos(c*cos(psai2));
    error = asin(c*sin(psai2))-acos(c*cos(psai2));
    Psai= [psai1; psai2];
    %
    %
    c1 = 1/r__(1)^5;
    c2 = 1/r__(2)^5;
    myFun2 = @(psai)[c1*sin(psai(1))+c2*sin(psai(2)) c1*cos(psai(1))+c2*cos(psai(2))];
    psai_0 = [30 01]*pi/180;
    [psai12,fval,exitflag,output] = fsolve(myFun2,psai_0);
%     [psai2,fval,exitflag,output] = fminsearch(@(psai) abs(myFun2(psai)),psai_0);
    error = myFun2(psai12);
    Psai = psai12;
end
%
function [Psai] = calcPsai3(x, y, MagPos, Psai1)

    for i=1:size(MagPos,1)
        x_ = MagPos(i,1);
        y_ = MagPos(i,2);
        z_ = MagPos(i,3);
        u_ = [x_ y_ z_];

        r_(i,:) = [x-x_ y-y_];
        r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
    end
    c32 = -r__(3)^5/r__(2)^5;
    c31 = -r__(3)^5/r__(1)^5;
    myFun = @(psai2)asin(c31*sin(Psai1)+c32*sin(psai2))-acos(c31*cos(Psai1)+c32*cos(psai2));
    psai2_0 = 30*pi/180;
    [psai2,fval,exitflag,output] = fsolve(myFun,psai2_0);
%     [psai2,fval,exitflag,output] = fminsearch(@(psai2) abs(myFun(psai2)),psai2_0);
    psai3 = asin(c31*sin(Psai1)+c32*sin(psai2));
    error = myFun(psai2);
    Psai = [Psai1; psai2; psai3];
    %
    %
    myFun2 = @(psai)asin(c31*sin(psai(1))+c32*sin(psai(2)))-acos(c31*cos(psai(1))+c32*cos(psai(2)));
    psai2_0 = [00; 00]*pi/180;
    [psai12,fval,exitflag,output] = fsolve(myFun2,psai2_0);
%     [psai12,fval,exitflag,output] = fminsearch(@(psai) abs(myFun2(psai)),psai2_0);
    psai3 = asin(c31*sin(psai12(1))+c32*sin(psai12(2)));
    Psai_ = [psai12(1); psai12(2); psai3];
    %
    %
    c1 = 1/r__(1)^5;
    c2 = 1/r__(2)^5;
    c3 = 1/r__(3)^5;
    myFun3 = @(psai)[c1*sin(psai(1))+c2*sin(psai(2))+c3*sin(psai(3)) c1*cos(psai(1))+c2*cos(psai(2))+c3*cos(psai(3))];
    psai_0 = [30 01 181]*pi/180;
    [psai123,fval,exitflag,output] = fsolve(myFun3,psai_0);
%     [psai2,fval,exitflag,output] = fminsearch(@(psai) abs(myFun3(psai)),psai_0);
    error = myFun3(psai123);
    Psai = psai123;
    %
    %
    %
%     figure
%     hold on
%     xx = min(MagPos(:,1)):0.01:max(MagPos(:,1));
%     yy = min(MagPos(:,2)):0.01:max(MagPos(:,2));
%     for k1=1:length(xx)
%         for k2=1:length(yy)
%            for i=1:size(MagPos,1)
%                 x_ = MagPos(i,1);
%                 y_ = MagPos(i,2);
%                 z_ = MagPos(i,3);            
% 
%                 r_(i,:) = [xx(k1)-x_ yy(k2)-y_];
%                 r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
%             end
%             counter = 0;
%             c1 = 1/r__(1,1)^5;
%             c2 = 1/r__(2,1)^5;
%             c3 = 1/r__(3,1)^5;
%             if c1+c2>c3
%                 counter = counter + 1;
%             end
%             if c1+c3>c2
%                 counter = counter + 1;
%             end
%             if c3+c2>c1
%                 counter = counter + 1;
%             end
%             if counter == 3
%                 plot(xx(k1), yy(k2), 'ro')
%             else
%                 plot(xx(k1), yy(k2), 'bo')
%             end
%         end
%     end
    
end
