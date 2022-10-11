function [eq_x,eq_y] = findEqPoints_Dynamics(x, y, Psai, maxTime) % Frho ro begire va ba delF hesab kone na inke dynamics bere psiCont ro brgire
    
    doPlot = 0;

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
            if ans1(end,j) == ans1(end-1,j)
                plot(ans1(end,j),ans1(end,n+j),'g.','MarkerSize', 16);
            end
        end
    end
    counter = 0;
    for j=1:n
        if ans1(end,j) == ans1(end-1,j)
            counter = counter + 1;
            eq_x(counter) = ans1(end,j); 
            eq_y(counter) = ans1(end,n+j);
        end
    end
    %
    %
    [C, ia, ic] = unique(round(eq_x,5));
    eq_x = eq_x(ia);
    eq_y = eq_y(ia);
    %
end