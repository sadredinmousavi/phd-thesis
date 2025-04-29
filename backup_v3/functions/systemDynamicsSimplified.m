function dydt = systemDynamicsSimplified(t,y,Psai)

n = size(y,1)/4;
m = 0;
[ind1, ind2, ind3, ind4, ind5] = makeIndiceForMRs(n);

%
x_mr = y(ind1);
y_mr = y(ind2);
xd_mr = y(ind3);
yd_mr = y(ind4);
%
fx = zeros(n+m,1);
fy = zeros(n+m,1);
dydt = zeros(size(y,1),1);

%
% Psai = psaiController(t);
mass = 4/3*pi*0.025^3*(3) / 1e3; %Rho_pla = 1.25 g/cm3 Rho_Neodymium = 7.5 g/cm3 r=0.025 cm
Cd = 16/3*0.01*(0.25/1e3);

% magnetic forces
for i=1:n
    F = force_field_symbolic(x_mr(i), y_mr(i), Psai);
    fx(i) = F(1);
    fy(i) = F(2);
end
% drag forces of MRs
for i=1:n
    fx(i) = fx(i) - Cd * xd_mr(i);
    fy(i) = fy(i) - Cd * yd_mr(i);
end

dydt(ind1) = y(ind3);
dydt(ind2) = y(ind4);
dydt(ind3) = fx(ind5,1)/mass;
dydt(ind4) = fy(ind5,1)/mass;




function [ind1, ind2, ind3, ind4, ind5] = makeIndiceForMRs(n)
    ind1 = 1:n;
    ind2 = n+1:2*n;
    ind3 = 2*n+1:3*n;
    ind4 = 3*n+1:4*n;
    ind5 = 1:n;
end

end