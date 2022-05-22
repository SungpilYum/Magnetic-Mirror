clear all
close all
m = 9.1e-31; %[kg]
I = 10e+11; %[current]
R0 = 1; %[m]
r = linspace(-10,10,100);
z = linspace(-20,10,100)';
[A, Br, Bz] = B_field_coil(r, z, I, R0);
f = figure;
filename = 'magnetic_mirror_result.gif';
f.Position = [100 100 700 800];
subplot(2,2,[1,3]);
contour(r, z, A, 200)
hold on
quiver(r, z, Br, Bz)
xlabel('x')
ylabel('z')
x = [1e-30,1e-30,-5]; %[x, y, z]
[A0, Br0, Bz0] = B_field_coil(0, -5, I, R0);
[Am, Brm, Bzm] = B_field_coil(0, -10, I, R0);
Rm = Bzm / Bz0;
Rm=Rm/10;
v0 = 10e+2;
v = [0, sqrt(v0.^2 / Rm), sqrt((1-1/Rm)*v0.^2)]; %[x, y, z]
Rm=Rm/0.001;
x2 = [-1e-30,-1e-30,-5]; %[x, y, z]
v2 = [0, sqrt(v0.^2 / Rm), sqrt((1-1/Rm)*v0.^2)]; %[x, y, z]
maxi = 30000;
for k = 1:maxi
    xx(k) = x(1);
    yy(k) = x(2);
    zz(k) = x(3);
    vx(k) = v(1);
    vy(k) = v(2);
    vz(k) = v(3);
    xx2(k) = x2(1);
    yy2(k) = x2(2);
    zz2(k) = x2(3);
    vx2(k) = v2(1);
    vy2(k) = v2(2);
    vz2(k) = v2(3);
    [A, Br, Bz] = B_field_coil(sqrt(x(1).^2 + x(2).^2), x(3), I, R0);
    B = [Br*(x(1)/sqrt(x(1).^2+x(2).^2)), Br*(x(2)/sqrt(x(1).^2+x(2).^2)), Bz];
    v_para_hist(k) = norm(B ./ norm(B) .* dot(v, B ./ norm(B)));
    v_perp_hist(k) = norm(v - B ./ norm(B) .* dot(v, B ./ norm(B)));
    mu_hist(k) = 1/2 * m * v_perp_hist(k).^2 / norm(B);
    E_hist(k) = 1/2 * m * norm(v).^2;
    [x, v] = pusher_boris(x, v, B);
    [A, Br2, Bz2] = B_field_coil(sqrt(x2(1).^2 + x2(2).^2), x2(3), I, R0);
    B2 = [Br2*(x2(1)/sqrt(x2(1).^2+x2(2).^2)), Br2*(x2(2)/sqrt(x2(1).^2+x2(2).^2)), Bz2];
    [x2, v2] = pusher_boris(x2, v2, B2);
    k;
end
subplot(2,2,[1,3]);
trajectory1 = plot(sqrt(xx(1).^2+yy(1).^2), zz(1), 'r.');
trajectory2 = plot(sqrt(xx2(1).^2+yy2(1).^2), zz2(1), 'b.');

for j = 1:400:maxi
    subplot(2,2,[1,3]);
    longss = 1000;
    if j > longss
        delete(trajectory1)
        trajectory1 = plot(sqrt(xx(j-longss:j).^2+yy(j-longss:j).^2), zz(j-longss:j), 'r.');
        delete(trajectory2)
        trajectory2 = plot(sqrt(xx2(j-longss:j).^2+yy2(j-longss:j).^2), zz2(j-longss:j), 'b.');
    else
        delete(trajectory1)
        trajectory1 = plot(sqrt(xx(1:j).^2+yy(1:j).^2), zz(1:j), 'r.');
        delete(trajectory2)
        trajectory2 = plot(sqrt(xx2(1:j).^2+yy2(1:j).^2), zz2(1:j), 'b.');
    end
    axis([-4 4 -16 6])
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
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
        if j == 1;
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
         imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    hold on;
end


hold off
