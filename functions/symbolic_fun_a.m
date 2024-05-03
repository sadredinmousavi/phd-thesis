function [] = symbolic_fun_a(values, fileName)
syms x y z m r mu_0 f real
syms m_norm m_norm_agent real

if nargin < 2
    fileName = 'dataTemp.mat';
end

%[rho theta phi m_norm mu_0]
% values = [
%     0.15 45*(pi/180) 000*(pi/180) 0.0127 4*pi*1e-7
%     0.15 45*(pi/180) 090*(pi/180) 0.0127 4*pi*1e-7
%     0.15 45*(pi/180) 180*(pi/180) 0.0127 4*pi*1e-7
%     0.15 45*(pi/180) 270*(pi/180) 0.0127 4*pi*1e-7
%     ];
for i=1:size(values,1)
    psai(i) = sym(strcat('psai', num2str(i)), 'real');
    array_{i} = num2str(i);
    str_psai{i} = sprintf('psai%d = psai(%d);\n', i, i);
    str_psai_sym{i} = sprintf('syms psai%d real;\n', i);
end

for i=1:size(values,1)
    rho(i) = sym(strcat('rho', num2str(i)), 'real');
end
for i=1:size(values,1)
    theta(i) = sym(strcat('theta', num2str(i)), 'real');
end
for i=1:size(values,1)
    phi(i) = sym(strcat('phi', num2str(i)), 'real');
end
%
waitText  = 'symbolic calculations - Please wait...';
waitIters = size(values,1);
waitHandle = waitbar(0,waitText);
for i=1:size(values,1)
    m_hat = [
        -cos( theta(i) )*cos( phi(i) )*cos( psai(i) )+sin( phi(i) )*sin( psai(i) )
        -cos( theta(i) )*sin( phi(i) )*cos( psai(i) )-cos( phi(i) )*sin( psai(i) )
        +sin( theta(i) )*cos( psai(i) )
        ];
    r_hat = [
        x - rho(i)*sin( theta(i) )*cos( phi(i) )
        y - rho(i)*sin( theta(i) )*sin( phi(i) )
        z - rho(i)*cos( theta(i) )
        ];
    r_norm = simplify(sqrt(r_hat(1)^2+r_hat(2)^2+r_hat(3)^2));
    r_hat = r_hat / simplify(sqrt(r_hat(1)^2+r_hat(2)^2+r_hat(3)^2));

%     simplify(m_hat(1)^2+m_hat(2)^2+m_hat(3)^2);
%     simplify(r_hat(1)^2+r_hat(2)^2+r_hat(3)^2);

    B = simplify(mu_0/4/pi*(m_norm/r_norm^3)*(3*(r_hat*r_hat')-eye(3))*m_hat);
    b_ = simplify( subs(B, [rho(i) theta(i) phi(i) m_norm mu_0], values(i,1:end-1)) );
    b_m_ = simplify( subs(m_norm_agent*B, [rho(i) theta(i) phi(i) m_norm mu_0 m_norm_agent], values(i,:)) );
    %
%     G = simplify(3*mu_0/4/pi*(m_norm/r_norm^4)*((m_hat*r_hat') + (r_hat*m_hat') - (5*(r_hat*r_hat')-eye(3))*(m_hat.*r_hat)));
%     ff = simplify(G(:,3));
%     f_ = simplify( subs(ff, [rho(i) theta(i) phi(i) m_norm mu_0], values(i,:)) );
%     % ff = simplify(G(2,3)-G(3,2));%ff = simplify(G(:,3)-G(3,:)');
%     % eval(subs(f_,[x y psai(1)], [05 05 1]))
%     % f__ = simplify( subs(f_, [x y], [0.05 0.05]) )
%     % fplot(f__(2), [-pi pi])
    if i<2
%         f = f_;
        b = b_;
        b_m = b_m_;
%         latex(f_(2))
    else
%         f = f + f_;
        b = b + b_;
        b_m = b_m + b_m_;
    end
    waitbar(i/waitIters,waitHandle, waitText);
end

cd('data')
save(fileName, 'b', 'b_m', 'str_psai', 'str_psai_sym')
cd('..')
makeFunctionsFromSymbolic(b, b_m, str_psai, str_psai_sym)

close(waitHandle)

end


