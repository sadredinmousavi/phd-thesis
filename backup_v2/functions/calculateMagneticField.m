function [B, Brho, Baxial] = calculateMagneticField(x,y,Psai)
    Npoints = length(x);
    Brho  = zeros(Npoints,Npoints); Baxial  = Brho;% Bz  = Brho;
    for a = 1:Npoints
        for b = 1:Npoints
            B_ = magnetic_field_symbolic(x(b), y(a), Psai);
            Brho(a,b)   = B_(1);
            Baxial(a,b) = B_(2);
            %Bz(a,b)     = B_(3);
        end
    end
    B = sqrt(Brho.^2+Baxial.^2);
end
