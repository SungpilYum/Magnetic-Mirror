function [x, v] = pusher_boris(x, v, B)

m = 9.1e-31;
q = 1.6e-19;
h = 4.3e-7;
E = [0, 0, 0];
% [A, Br, Bz] = B_field_coil(sqrt(x(1).^2 + x(2).^2), x(3));
% B = [Br*(x(1)/sqrt(x(1).^2+x(2).^2)), Br*(x(2)/sqrt(x(1).^2+x(2).^2)), Bz];

v_minus = double.empty(3,0);
v_prime = double.empty(3,0);
v_plus =  double.empty(3,0);
v_next =  double.empty(3,0);

v_minus(1) = v(1) + (q * h * 0.5 * E(1)) / m;
v_minus(2) = v(2) + (q * h * 0.5 * E(2)) / m;
v_minus(3) = v(3) + (q * h * 0.5 * E(3)) / m;

v_prime(1) = v_minus(1) + (v_minus(2) * B(3) - v_minus(3) * B(2)) * (q * h * 0.5) / m;
v_prime(2) = v_minus(2) + (v_minus(3) * B(1) - v_minus(1) * B(3)) * (q * h * 0.5) / m;
v_prime(3) = v_minus(3) + (v_minus(1) * B(2) - v_minus(2) * B(1)) * (q * h * 0.5) / m;

v_plus(1) = v_minus(1) + (v_prime(2) * B(3) - v_prime(3) * B(2)) * (4 * q * h * m) / (q * q * (B(1).^2 + B(2).^2 + B(3).^2) * h * h + 4 * m * m);
v_plus(2) = v_minus(2) + (v_prime(3) * B(1) - v_prime(1) * B(3)) * (4 * q * h * m) / (q * q * (B(1).^2 + B(2).^2 + B(3).^2) * h * h + 4 * m * m);
v_plus(3) = v_minus(3) + (v_prime(1) * B(2) - v_prime(2) * B(1)) * (4 * q * h * m) / (q * q * (B(1).^2 + B(2).^2 + B(3).^2) * h * h + 4 * m * m);

v_next(1) = v_plus(1) + (q * h * 0.5 * E(1)) / m;
v_next(2) = v_plus(2) + (q * h * 0.5 * E(2)) / m;
v_next(3) = v_plus(3) + (q * h * 0.5 * E(3)) / m;

x(1) = x(1) + v_next(1) * h;
x(2) = x(2) + v_next(2) * h;
x(3) = x(3) + v_next(3) * h;

v(1) = v_next(1);
v(2) = v_next(2);
v(3) = v_next(3);


end