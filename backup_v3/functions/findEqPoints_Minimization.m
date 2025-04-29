function [eq_x,eq_y] = findEqPoints_Minimization(x, y, Psai, stepSize) % Frho ro begire va ba delF hesab kone na inke dynamics bere psiCont ro brgire
    % steepest descent
    
    doPlot = 0;
    
    if nargin < 4
        stepSize = ( max(x)-min(x) ) / 30;
    end
    maxIter = 500;
    tol = 1e-8;
    
    stepSizeRegion = ( max(x)-min(x) ) / 5;
    x__ = min(x)+0.0:stepSizeRegion:max(x)-0.0;
    y__ = min(y)+0.0:stepSizeRegion:max(y)-0.0;
    [x_,y_] = meshgrid(x__,y__);
    x_ = reshape(x_, 1, []);
    y_ = reshape(y_, 1, []);
%     x_ = 0.1;
%     y_ = 0.1;
    %
% %     ans0 = [x_ y_ zeros(1,length(x_)) zeros(1,length(x_))];
% %     [t,ans1] = ode45(@(t,y) systemDynamicsSimplified(t,y,Psai), [0 stepSize], ans0);
% %     n = size(ans1,2)/4;
    %
    ans0 = [x_ y_];
    ans1 = ans0;
    n = length(x_);
    Frho = zeros(1, n); Faxial_ = Frho;
    [~, Frho_, Faxial_] = calculateForceField(x,y,Psai);
    alpha = stepSize * ones(1,n) / max(max(Frho_));
    for i=1:maxIter
        for j = 1:n
            if min(x) < ans1(i,j) &&  ans1(i,j) < max(x) && min(y) < ans1(i,n+j) &&  ans1(i,n+j) < max(y)
                F = force_field_symbolic(ans1(i,j), ans1(i,n+j), Psai);
            else
                F = [0;0];
            end
            Frho(1,j) = F(1);
            Faxial(1,j) = F(2);
        end
        ans1(i+1,:) = ans1(i,:) + [alpha.*Frho alpha.*Faxial];
    end
    %
    %
    %
    if doPlot
        [~, Frho, Faxial] = calculateForceField(x,y,Psai);
        plot(ans1(end,1:n),ans1(end,n+1:2*n),'b.','MarkerSize', 16);
        hold on
        for j=1:n
            plot(ans1(:,j),ans1(:,n+j),'r-', 'LineWidth',1 );
        end
        plot_field = streamslice(x,y,Frho,Faxial,'method','cubic');
        set(plot_field,'Color','black','LineWidth',1.2);
        %
        for j=1:n
            %if ans1(end,j) - ans1(end-1,j) < tol
            if ans1(end,j) == ans1(end-1,j)
                if min(x) < ans1(end,j) &&  ans1(end,j) < max(x) && min(y) < ans1(end,n+j) &&  ans1(end,n+j) < max(y)
                    plot(ans1(end,j),ans1(end,n+j),'g.','MarkerSize', 16);
                end
            end
        end
    end
    counter = 0;
    for j=1:n
        %if ans1(end,j) - ans1(end-1,j) < tol
        if ans1(end,j) == ans1(end-1,j)
            if min(x) < ans1(end,j) &&  ans1(end,j) < max(x) && min(y) < ans1(end,n+j) &&  ans1(end,n+j) < max(y)
                counter = counter + 1;
                eq_x(counter) = ans1(end,j); 
                eq_y(counter) = ans1(end,n+j);
            end
        end
    end
    %
    %
    try
        [C, ia, ic] = unique(round(eq_x,5));
        eq_x = eq_x(ia);
        eq_y = eq_y(ia);
    catch
        eq_x = 0;
        eq_y = 0;
    end
    %
end