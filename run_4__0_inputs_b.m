%% Cylinder magnetic field visualization
 % This script is a DEMO for the visualization of the magnetic flux density 
 % field lines of an axially magnetized cylinder. 

clc, close all, clearvars
addpath(genpath('functions'))
addpath(genpath('inputFiles'))


inputFile001_7


% [r,isInContact,normalVector] = calculateDistanceToWallLinear([0;0], [-3;0;0;3])

%% Define parameters
vars = args;
vars.args = args;


%% Define PMs locations
% values --> [x y z m_norm mu_0 m_agents]  for PMs
MagPos = [values(:,1:3) ones(size(values,1),1) zeros(size(values,1),2)];
for i=1:size(MagPos, 1)
    MagPos(i,4:6) = round(MagPos(i,4:6) ./ norm(MagPos(i,4:6)), 5);
end
%
%
vars.values = values;
vars.MagPos = MagPos;
vars.phi = phi;
%
vars.calcPsaiFromEqFunc = calcPsaiFromEqFunc;
vars.findEqFromMinimization = findEqFromMinimization;


%% Define filing configs
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
    end
end
vars.fileName1 = fileName1;
vars.fileName2 = fileName2;
vars.fileName3 = fileName3;
cd('..')
%
%
%
dlgTitle    = 'PreCalculations';
dlgQuestion = 'Do you want symbolic calculations to be done ?';
choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
if strcmp(choice, 'Yes') == 1
    DoPreCalc = 1;
    symbolicFunctionHandle(values, fileName3, vars);
else
    DoPreCalc = 0;
end
vars.DoPreCalc = DoPreCalc;
% makeFunctionsFromSymbolic(b, b_m, str_psai, str_psai_sym) % inside _datas


%% Define space and plot parameters and options
vars.x_space = x_space;
vars.y_space = y_space;
vars.z_space = z_space;
vars.plotDomain = plotDomain;
%
vars.plotOptions.dynamic = plotOptionsDyn;


%% Define MRs locations
vars.x_mr_0 = x_mr_0;
vars.y_mr_0 = y_mr_0;
vars.r_mr_0 = r_mr_0;

%% Define FPs locations
if exist('x_fp_0', 'var') == 1
    vars.x_fp_0 = x_fp_0;
    vars.y_fp_0 = y_fp_0;
    vars.t_fp_0 = t_fp_0;
    vars.fps = fps;
else
    vars.x_fp_0 = [];
    vars.y_fp_0 = [];
    vars.r_fp_0 = [];
    vars.fps = [];
end

%% Define Walls
if exist('walls', 'var') == 1
    vars.walls = walls;
else
    vars.walls = [];
end


%% calculate Psai_0 in order to have equilibrium in 2 points (using 4 or 6 magnets)
dlgTitle    = 'Hand Calculations';
dlgQuestion = 'Do you want eqPoint calculations to be done ?';
choice = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes');
if strcmp(choice, 'Yes') == 1
    % eqPoint1 = [ [0;+0.05;+0.00] [300;+0.03;+0.04] [600;+0.03;+0.09] ];
    % eqPoint2 = [ [0;-0.11;-0.02] [300;-0.10;-0.09] [600;+0.07;-0.10] ];
    eq_points{1} = [+0.05 +0.03];
    eq_points{2} = [+0.05 +0.03];
%     eq_points{3} = [+0.00 -0.08];
    vars.plotOptions.static.plotEqPoints = 1;
    vars.plotOptions.static.printLambdaValues = 1;
    vars.plotOptions.static.printPsaiValues = 1;
    vars.plotOptions.static.eq_points = eq_points;
    

    [rankM, error, hasAns, isStable, Psai, hessian, otherOutputs] = vars.calcPsaiFromEqFunc(eq_points, MagPos);
    printFig2(vars, Psai, 0);
%     PsaiSerie = otherOutputs.PsaiSerie;
%     for cnt=1:length(PsaiSerie)
%         Psai = PsaiSerie(:,cnt);
%         printFig(vars, Psai, 0);
%     end
    
    %
    [Fx1, Fy1, Fz1] = CylFfield3(Psai,eq_points{1});
    [Bx1, By1, Bz1] = CylBfield3(Psai,eq_points{1});
    %
    [Fx2, Fy2, Fz2] = CylFfield3(Psai,eq_points{2});
    [Bx2, By2, Bz2] = CylBfield3(Psai,eq_points{2});
    %
%     [Fx3, Fy3, Fz3] = CylFfield3(Psai,eq_points{3});
%     [Bx3, By3, Bz3] = CylBfield3(Psai,eq_points{3});
    %
    printFig(vars, Psai, 0);
    %
%     [V,D] = eig(hessian.point1); % V(:,i)
%     angle1_point1 = atan(V(2,1)/V(1,1))*(180/pi);
%     angle2_point1 = atan(V(2,2)/V(1,2))*(180/pi);
%     d1_point1 = D(1,1);
%     d2_point1 = D(2,2);
%     [V,D] = eig(hessian.point2); % V(:,i)
%     angle1_point2 = atan(V(2,1)/V(1,1))*(180/pi);
%     angle2_point2 = atan(V(2,2)/V(1,2))*(180/pi);
%     d1_point2 = D(1,1);
%     d2_point2 = D(2,2);



    newStep = 0.01;
    figure
    hold on
    plot(MagPos(:,1), MagPos(:,2), 'r.', 'MarkerSize', 15)
    plot(eq_point1(1), eq_point1(2), '+r')
    for x=-domain:newStep:domain
        for y=-domain:newStep:domain
            eq_points{1} = eq_point1;
            eq_points{2} = [x y];
            [rankM, error, hasAns, isStable, Psai, hessian] = vars.calcPsaiFromEqFunc(eq_points, MagPos);
    %         if hasAns
            if isStable && hasAns
                plot(x, y, 'co')
            else
                plot(x, y, 'ko')
            end
            if hasAns
                a=1;
            end
        end
    end
    axis square
end


%% calculate Psai_0 in order to have equilibrium in 3 points (using 6 magnets)
% eq_points{1} = [-0.02 +0.04];%eq_point1 = [-0.01 +0.05];
% eq_points{2} = [+0.04 -0.11];%eq_point2 = [+0.04 -0.13];
% eq_points{3} = [+0.04 -0.11];%eq_point3 = [+0.06 -0.11];
% 
% [rankM, error, hasAns, isStable, Psai] = check_6PM3(eq_points, MagPos);
% printFig(vars, Psai, 0);
% newStep = 0.01;
% 
% allPointCount = 0;
% equalPointCount = 0;
% for x2=-domain:newStep:domain
%     for y2=-domain:newStep:domain
%         x3=x2;y3=y2;
% %         for x3=-domain:newStep:domain
% %             for y3=-domain:newStep:domain
%                 [rankM,error, hasAns, isStable, Psai] = check_6PM3(eq_point1, [x2 y2], [x3 y3], MagPos);
% %                 if hasAns
%                 if isStable && hasAns
%                     allPointCount = allPointCount + 1;
%                     eqpoint1(allPointCount, :) = eq_point1;
%                     eqpoint2(allPointCount, :) = [x2 y2];
%                     eqpoint3(allPointCount, :) = [x3 y3];
%                     if x2==x3 && y2==y3
%                         equalPointCount = equalPointCount + 1;
%                     end
%                 end
%                 if hasAns
%                     a=1;
%                 end
% %             end
% %         end
%     end
% end
% %
% figure
% hold on
% plot(MagPos(:,1), MagPos(:,2), 'r.', 'MarkerSize', 15)
% plot(eq_point1(1), eq_point1(2), '+k')
% for x=-domain:newStep:domain
%     for y=-domain:newStep:domain
%         plot(x, y, 'ko')
%     end
% end
% for i=1:length(eqpoint2)
%     plot(eqpoint2(i,1), eqpoint2(i,2), 'co')
% end


%% Define path and eqPoints desired locations during simulation time
% path inputs --> [point1 point2 startTime endTime dt]
% % eqPoint = [t1 t2 t3 ...; x1 x2 x3 ...; y1 y2 y3 ...]
vars.tspan = tspan;
vars.eqPoints = eqPoints;
vars.usePrepaidPsai = usePrepaidPsai;
if usePrepaidPsai
    vars.Psai = Psai;
end
% % [eq_x,eq_y] = findEqPoints_Dynamics(x_space,y_space,Psai_0);
% [eq_x,eq_y] = findEqPoints_Minimization(x_space,y_space,Psai_0);
% vars.eqPoints.x_eqPoints_0 = eq_x;
% vars.eqPoints.y_eqPoints_0 = eq_y;

designPsaiController(vars);


%%
cd('data')
save(fileName1, '-struct', 'vars')
cd('..')

run_4__1



%%
function [eqPoint1] = defineLinearPath(eqPoint0, eqPoint1, t0, refreshTime, speed)
    point1_x0 = eqPoint0(1);
    point1_y0 = eqPoint0(2);
    point1_x1 = eqPoint1(1);
    point1_y1 = eqPoint1(2);
    %
    speed_x = speed * (point1_x1-point1_x0)/sqrt((point1_x1-point1_x0)^2+(point1_y1-point1_y0)^2);
    speed_y = speed * (point1_y1-point1_y0)/sqrt((point1_x1-point1_x0)^2+(point1_y1-point1_y0)^2);
    t1 = max( abs( (point1_x1-point1_x0) / speed_x), abs( (point1_y1-point1_y0) / speed_y) );
    time = 0:refreshTime:t1;
    eqPoint1 = zeros(3,length(time));
    for i=1:length(time)
        eqPoint1(:,i) = [t0+time(i); point1_x0+speed_x*time(i); point1_y0+speed_y*time(i)];
    end
end
%
function [eqPoint1] = defineCirclePath(r, deg, t0, refreshTime, speed)
    pointFunCyln = @(t,r,theta)[t; r*cos(theta*(pi/180)); r*sin(theta*(pi/180))];
    omega = speed/r*(180/pi); %deg/s
    t1 = abs( (deg(2)-deg(1)) / omega );
    time = 0:refreshTime:t1;
    eqPoint1 = zeros(3,length(time));
    for i=1:length(time)
        eqPoint1(:,i) = pointFunCyln(t0+time(i),r,deg(1)+omega*time(i));
    end
end
%
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
    %
    %
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

    
end
%
function [Psai_0,values] = calcPsaiNew(x, y, values, a, phi, z)
    input_0 = [a 0.1 0.1 0.1];
    input_0 = [a(1) 0.1 0.1 0.1];
    opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
    [ans,fval,exitflag,output] = fsolve(@(input)myFun0(input, phi, x, y, z),input_0,opts);
    a = ans(1:3)';
    Psai_0 = ans(4:6)';
    values_ = [
        a(1)*cos(phi(1)) a(1)*sin(phi(1))
        a(2)*cos(phi(2)) a(2)*sin(phi(2))
        a(3)*cos(phi(3)) a(3)*sin(phi(3))
    ];
    values = [values_ values(:,3:end)];
end
%
function [ans] = myFun0(input, phi, x, y, z)
%     a(1) = input(1);
%     a(2) = input(2);
%     a(3) = input(3);
    a = [input(1) 0.25 0.25];
    psai(1) = input(2);
    psai(2) = input(3);
    psai(3) = input(4);
    for i=1:3
        x_ = a(i)*cos(phi(i));
        y_ = a(i)*sin(phi(i));
        z_ = z;

        r_(i,:) = [x-x_ y-y_];
        r__(i,1) = sqrt( r_(i,1)^2 + r_(i,2)^2 + (0-z_)^2 );
    end
    c1 = 1/r__(1)^5;
    c2 = 1/r__(2)^5;
    c3 = 1/r__(3)^5;
    ans = [
        c1*sin(psai(1))+c2*sin(psai(2))+c3*sin(psai(3))
        c1*cos(psai(1))+c2*cos(psai(2))+c3*cos(psai(3))
        ];
end
