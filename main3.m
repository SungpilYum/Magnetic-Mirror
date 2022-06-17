clear all
close all
m = 9.1e-31; %[kg]
I = 10e+11; %[current]
r0 = 2; %[m], minor radius
R0 = 3; %[m], major radius
x = linspace(-10,10,200);
y = linspace(-10,10,200);
z = linspace(-10,10,200);
[X, Y, Z] = meshgrid(x, y, z);
N = 16;%number of toroidal coils
[A, Bx, By, Bz] = B_field_torus(X, Y, Z, I, r0, R0, N);
[sx, sy, sz] = meshgrid(linspace(0,8,3), linspace(0,8,3), linspace(-3,3,1));
% plot_B_field = streamline(X,Y,Z,Bx,By,Bz, sx, sy, sz);
% view(90, 0)
figure(3)
plot_B_field2 = streamslice(X,Y,Z,Bx,By,Bz, [], [], 0, 5,"arrows","cubic");
view(3)
for i = 1:N
    phi = i*2*pi/N;
    hold on
    plotCircle3D([R0*cos(phi),R0*sin(phi),0], [-R0*sin(phi),R0*cos(phi),0], r0);
end
hold off
pause

xlabel('x')
ylabel('y')
ylabel('z')

x = [R0+r0/2,1e-30,1e-30]; %[x, y, z]
[Am, Bxm, Bym, Bzm] = B_field_torus(R0-r0, 0, 0, I, r0, R0, N);
[A0, Bx0, By0, Bz0] = B_field_torus(R0+r0, 0, 0, I, r0, R0, N);
Rm = sqrt(Bxm.^2 + Bym.^2 + Bzm.^2) / sqrt(Bx0.^2 + By0.^2 + Bz0.^2);
v0 = 10e+2;
Rm = Rm / 10;
v = [0, sqrt(v0.^2 / Rm), sqrt((1-1/Rm)*v0.^2)]; %[x, y, z]
maxi = 50000;
for k = 1:maxi
    xx(k) = x(1);
    yy(k) = x(2);
    zz(k) = x(3);
    vx(k) = v(1);
    vy(k) = v(2);
    vz(k) = v(3);
    [A, Bx, By, Bz] = B_field_torus(x(1), x(2), x(3), I, r0, R0, N);
    B = [Br*(x(1)/sqrt(x(1).^2+x(2).^2)), Br*(x(2)/sqrt(x(1).^2+x(2).^2)), Bz];
    v_para_hist(k) = norm(B ./ norm(B) .* dot(v, B ./ norm(B)));
    v_perp_hist(k) = norm(v - B ./ norm(B) .* dot(v, B ./ norm(B)));
    mu_hist(k) = 1/2 * m * v_perp_hist(k).^2 / norm(B);
    E_hist(k) = 1/2 * m * norm(v).^2;
    [x, v] = pusher_boris(x, v, B);
    [A, Br2, Bz2] = B_field_coil(sqrt(x2(1).^2 + x2(2).^2), x2(3), I, r0);
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
