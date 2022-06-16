function [Bx, By, Bz] = B_field_toros(x, y, z, I, r0, R0)
nu = 8.85418782e-12;
A_  = @(a, I, r, z) ...
    (nu./(4.*pi)) .* ((pi.*a.^2 .* I .* r)./(a.^2+r.^2+z.^2).^(3./2)) .* (1+((15.*a.^2.*r.^2)./(8.*(a.^2+r.^2+z.^2).^2)));
Br_ = @(a, I, r, z) ...
    (15.*I.*a.^4.*nu.*r.^3.*z)./(8.*(a.^2 + r.^2 + z.^2).^(9./2)) + (3.*I.*a.^2.*nu.*r.*z.*((15.*a.^2.*r.^2)./(8.*(a.^2 + r.^2 + z.^2).^2) + 1))./(4.*(a.^2 + r.^2 + z.^2).^(5./2));
Bz_ = @(a, I, r, z) ...
    (I.*a.^2.*nu.*((15.*a.^2.*r.^2)./(8.*(a.^2 + r.^2 + z.^2).^2) + 1))./(4.*(a.^2 + r.^2 + z.^2).^(3./2)) - (3.*I.*a.^2.*nu.*r.^2.*((15.*a.^2.*r.^2)./(8.*(a.^2 + r.^2 + z.^2).^2) + 1))./(4.*(a.^2 + r.^2 + z.^2).^(5./2)) - (I.*a.^2.*nu.*r.*((15.*a.^2.*r.^3)./(2.*(a.^2 + r.^2 + z.^2).^3) - (15.*a.^2.*r)./(4.*(a.^2 + r.^2 + z.^2).^2)))./(4.*(a.^2 + r.^2 + z.^2).^(3./2));


Bx = zeros(size(x));
By = zeros(size(y));
Bz = zeros(size(z));

for i = 1:8
phi = i * 2*pi / 8;
a = r0;
% Z = y;
% R = sqrt((x-R0).^2 + z.^2);
% Br0 = Br_(a, I, R, Z);
% Bx1 = (x-R0).*sqrt(1./((x-R0).^2 + z.^2)).*Br0;
% By1 = Bz_(a, I, R, Z);
% Bz1 = z.*sqrt(1./((x-R0).^2 + z.^2)).*Br0;

Z = y.*cos(phi) - x.*sin(phi);
R = sqrt((y.*sin(phi) + x.*cos(phi) - R0).^2 + z.^2);
Br0 = Br_(a, I, R, Z);
Bz0 = Bz_(a, I, R, Z);
Bx1 = Bz0.* -sin(phi) + cos(phi).*(y.*sin(phi) + x.*cos(phi) - R0).*sqrt(1./((y.*sin(phi) + x.*cos(phi) - R0).^2 + z.^2)).*Br0;
By1 = Bz0.*  cos(phi) + sin(phi).*(y.*sin(phi) + x.*cos(phi) - R0).*sqrt(1./((y.*sin(phi) + x.*cos(phi) - R0).^2 + z.^2)).*Br0;
Bz1 = z.*sqrt(1./((y.*sin(phi) + x.*cos(phi) - R0).^2 + z.^2)).*Br0;

Bx = Bx + Bx1;
By = By + By1;
Bz = Bz + Bz1;
end