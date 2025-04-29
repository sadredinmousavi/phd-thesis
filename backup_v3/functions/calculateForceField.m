function [F, Frho, Faxial] = calculateForceField(x,y,Psai)
    Npoints = length(x);
    Frho  = zeros(Npoints,Npoints); Faxial  = Frho;% Fz  = Frho;
    for a = 1:Npoints
        for b = 1:Npoints
            F_ = force_field_symbolic(x(b), y(a), Psai);
            Frho(a,b)   = F_(1);
            Faxial(a,b) = F_(2);
            %Fz(a,b)     = F_(3);
        end
    end
    F = sqrt(Frho.^2+Faxial.^2);
end

