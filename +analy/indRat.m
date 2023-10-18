function [ratio] = indRat(f,a)
%INDRAT 
%   a is wire radius.

sigma = 5.8e7; %conductivitiy of the wire
mu_0 = 4*pi*1e-7;
w = 2 * pi *f;
k = (sigma * w * mu_0).^0.5;
[r0 i0] = be(0,k*a);
[r1 i1] = bep(k*a);
[r2 i2] = be(2,k*a);

ratio = -2*pi*a/sigma .*k .* (r2.*r1+i2.*i1)./(r0.^2+i0.^2);

    function [ber bei] = be(nu,z)
        val = besselj(nu,z.*exp(3i*pi/4));
        ber = real(val);
        bei = imag(val);
    end

    function [ber bei] = bep(z)
        [ber1, bei1] = be(1,z);
        ber = (ber1+bei1)./sqrt(2);
        bei = (-ber1+bei1)./sqrt(2);
    end
end

