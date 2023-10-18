function [ratio] = skinRat(f,a)
%SKINRAT 
%   a is wire radius

sigma = 5.8e7; %conductivitiy of the wire
mu_0 = 4*pi*1e-7;

k = (sigma * 2*pi*f * mu_0).^0.5;
[m0 t0] = M_0(k*a);
[m1 t1] = M_1(k*a);

ratio = k.*a/2 .* m0 ./ m1 .* cos(t0 - t1 + 0.75*pi);

    function [M0 theta_0] = M_0(z)
        J0 = besselj(0,z*(1i^1.5));
        [theta_0 M0] = cart2pol(real(J0),imag(J0));
    end
    
    function [M1 theta_1] = M_1(z)
        J1 = besselj(1,z*(1i^1.5));
        [theta_1 M1] = cart2pol(real(J1),imag(J1));
    end
end

