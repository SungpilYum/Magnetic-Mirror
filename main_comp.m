clear all
close all
m = 9.1e-31; %[kg]
I = 10e+5; %[current]
R0 = 1; %[m]
Rc = 0.7; %[m]
r = linspace(-10,10,100);
z = linspace(-40,10,100)';
[A, Br, Bz] = B_field_compressor(r, z, I, R0, Rc);
f = figure;
f.Position = [100 100 700 800];
subplot(2,2,[1,3]);
contour(r, z, A, 200)
hold on
quiver(r, z, Br, Bz)
xlabel('x')
ylabel('z')
x = [1e-30,1e-30,-5]; %[x, y, z]
[A0, Br0, Bz0] = B_field_compressor(0, -5, I, R0, Rc);
[Am, Brm, Bzm] = B_field_compressor(0, -10, I, R0, Rc);
Rm = Bzm / Bz0;
Rm = Rm/0.01
v0 = 10e+2;
v = [0, sqrt(v0.^2 / Rm), -sqrt((1-1/Rm)*v0.^2)]; %[x, y, z]
maxi = 5000;
for k = 1:maxi
    xx(k) = x(1);
    yy(k) = x(2);
    zz(k) = x(3);
    vx(k) = v(1);
    vy(k) = v(2);
    vz(k) = v(3);
    I = I*(1+1*sin(k/100*pi));
    [A, Br, Bz] = B_field_compressor(sqrt(x(1).^2 + x(2).^2), x(3), I, R0, Rc);
    B = [Br*(x(1)/sqrt(x(1).^2+x(2).^2)), Br*(x(2)/sqrt(x(1).^2+x(2).^2)), Bz];
    v_para_hist(k) = norm(B ./ norm(B) .* dot(v, B ./ norm(B)));
    v_perp_hist(k) = norm(v - B ./ norm(B) .* dot(v, B ./ norm(B)));
    mu_hist(k) = 1/2 * m * v_perp_hist(k).^2 / norm(B);
    E_hist(k) = 1/2 * m * norm(v).^2;
    [x, v] = pusher_boris(x, v, B);
    k;
end
subplot(2,2,[1,3]);
trajectory1 = plot(sqrt(xx(1).^2+yy(1).^2), zz(1), 'r.');
hold on
for j = 1:20:maxi
    subplot(2,2,[1,3]);
    I = I*(1+1*sin(j/100*pi));
    [A, Br, Bz] = B_field_compressor(r, z, I, R0, Rc);
    contour(r, z, A, 200)
    hold on
    quiver(r, z, Br, Bz)
    pause
    if j > 100
        delete(trajectory1)
        trajectory1 = plot(sqrt(xx(j-100:j).^2+yy(j-100:j).^2), zz(j-100:j), 'r.');
    else
        delete(trajectory1)
        trajectory1 = plot(sqrt(xx(1:j).^2+yy(1:j).^2), zz(1:j), 'r.');
    end
    axis([-5.5 5.5 -30 6])
    pause(0.1)
    subplot(2,2,2)
    hold on
    grid on
    plot(j, v_para_hist(j), 'g.', j, v_perp_hist(j), 'b.')
    xlim([0 maxi])
    ylim([-10e+1 v0+10e+1])
    legend('|v_{||}|', '|v_\perp|');
    hold off

    subplot(2,2,4)
    grid on
    hold on
    plot(j, mu_hist(j),'r.', j, E_hist(j), 'mo', j, v_para_hist(j), 'g.', j, v_perp_hist(j), 'b.')
    xlim([0 maxi])
    ylim([-1e-10 1e-10])
    legend('\mu', 'Energy');
    hold off
end


hold off
