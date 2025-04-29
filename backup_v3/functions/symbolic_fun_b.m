function [output] = symbolic_fun_b(epms, args)
syms x y z m r mu_0 f real
syms m_norm m_norm_agent real

if nargin < 2
    fileName = 'dataTemp.mat';
end

%[xx yy zz m_norm mu_0]
for i=1:epms.magNum
    psai(i) = sym(strcat('psai', num2str(i)), 'real');
    array_{i} = num2str(i);
    str_psai{i} = sprintf('psai%d = psai(%d);\n', i, i);
    str_psai_sym{i} = sprintf('syms psai%d real;\n', i);
end

for i=1:epms.magNum
    xx(i) = sym(strcat('xx', num2str(i)), 'real');
end
for i=1:epms.magNum
    yy(i) = sym(strcat('yy', num2str(i)), 'real');
end
for i=1:epms.magNum
    zz(i) = sym(strcat('zz', num2str(i)), 'real');
end
%
waitText  = 'symbolic calculations - Please wait...';
waitIters = epms.magNum;
waitHandle = waitbar(0,waitText);
for i=1:epms.magNum
    m_hat = [
        - sin( psai(i) ) * sin( epms.phi(i) )
        sin( psai(i) ) * cos( epms.phi(i) )
        cos( psai(i) )
        ];
    r_hat = [
        x - xx(i)
        y - yy(i)
        z - zz(i)
        ];
    r_norm = simplify(sqrt(r_hat(1)^2+r_hat(2)^2+r_hat(3)^2));
    r_hat = r_hat / r_norm;

%     simplify(m_hat(1)^2+m_hat(2)^2+m_hat(3)^2);
%     simplify(r_hat(1)^2+r_hat(2)^2+r_hat(3)^2);

    B = simplify(mu_0/4/pi*(m_norm/r_norm^3)*(3*(r_hat*r_hat')-eye(3))*m_hat);
    b_ = simplify( subs(B, [xx(i) yy(i) zz(i) m_norm mu_0], [epms.MagPos(i,1:3) args.pm.m args.mu_0]) );
    b_m_ = simplify( subs(m_norm_agent*B, [xx(i) yy(i) zz(i) m_norm mu_0 m_norm_agent], [epms.MagPos(i,1:3) args.pm.m args.mu_0 args.mr.m]) );
    f_ = 3*mu_0/4/pi*(m_norm*m_norm_agent/r_norm^4)*cos(psai(i))*r_hat(1:2,1);
    f_newformula_ = simplify( subs(f_, [xx(i) yy(i) zz(i) m_norm mu_0 m_norm_agent], [epms.MagPos(i,1:3) args.pm.m args.mu_0 args.mr.m]) );
    %
%     G = simplify(3*mu_0/4/pi*(m_norm/r_norm^4)*((m_hat*r_hat') + (r_hat*m_hat') - (5*(r_hat*r_hat')-eye(3))*(m_hat.*r_hat)));
%     g_ = simplify( subs(G, [xx(i) yy(i) zz(i) m_norm mu_0], [epms.MagPos(i,1:3) args.pm.m args.mu_0]) );
%     g_m_ = simplify( subs(m_norm_agent*G, [xx(i) yy(i) zz(i) m_norm mu_0 m_norm_agent], [epms.MagPos(i,1:3) args.pm.m args.mu_0 args.mr.m]) );
% %     G_3 = simplify(G(:,3));
% %     g_3 = simplify( subs(G_3, [rho(i) theta(i) phi(i) m_norm mu_0], values(i,:)) );
% %     % ff = simplify(G(2,3)-G(3,2));%ff = simplify(G(:,3)-G(3,:)');
% %     % eval(subs(f_,[x y psai(1)], [05 05 1]))
% %     % f__ = simplify( subs(f_, [x y], [0.05 0.05]) )
% %     % fplot(f__(2), [-pi pi])
    if i<2
        f_newformula = f_newformula_;
        b = b_;
        b_m = b_m_;
%         g = g_;
%         g_m = g_m_;
% %         latex(f_(2))
    else
        f_newformula = f_newformula + f_newformula_;
        b = b + b_;
        b_m = b_m + b_m_;
%         g = g + g_;
%         g_m = g_m + g_m_;
    end
    waitbar(i/waitIters,waitHandle, waitText);
end

output.b = b;
output.b_m = b_m;
output.f = f_newformula;
% output.g = g;
% output.g_m = g_m;
output.str_psai = str_psai;
output.str_psai_sym = str_psai_sym;

% makeFunctionsFromSymbolic(b, b_m, str_psai, str_psai_sym)

close(waitHandle)

end


