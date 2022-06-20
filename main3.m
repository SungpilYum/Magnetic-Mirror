clear all
close all
% m = 9.1e-31; %[kg] electron
m = 6.644657230e-27; %[kg] alpha particle
I = 50e+12; %[current]
% r0 = 1; %[m], minor radius
% R0 = 1.8; %[m], major radius
r0 = 4; %[m], minor radius
R0 = 15; %[m], major radius
x = linspace(-R0-r0, R0+r0, 100);
y = linspace(-R0-r0, R0+r0, 100);
z = linspace(-r0, r0, 100);
[X, Y, Z] = meshgrid(x, y, z);

N = 16 ;  %number of toroidal coils

[A, Bx, By, Bz] = B_field_torus(X, Y, Z, I, r0, R0, N);
% [sx, sy, sz] = meshgrid(linspace(0,R0+r0,3), linspace(0,R0+r0,3), linspace(-r0,r0,1));
figure(1)
subplot(2,2,[1,2]);
hold on
plot_B_field2 = streamslice(X,Y,Z,Bx,By,Bz, [], [], 0, 1,"arrows","cubic");
view(3)
for i = 1:N
    phi = i*2*pi/N;
    plotCircle3D([R0*cos(phi),R0*sin(phi),0], [-R0*sin(phi),R0*cos(phi),0], r0);
end
xlabel('x')
ylabel('y')
zlabel('z')
x = [R0+r0/3,1e-30,1e-30]; %[x, y, z]
[Am, Bxm, Bym, Bzm] = B_field_torus(R0-r0, 0, 0, I, r0, R0, N);
[A0, Bx0, By0, Bz0] = B_field_torus(R0+r0, 0, 0, I, r0, R0, N);
Rm = sqrt(Bxm.^2 + Bym.^2 + Bzm.^2) / sqrt(Bx0.^2 + By0.^2 + Bz0.^2);
% v0 = 50e+5;
% K = 0.5; % electron energy [eV]
K =1e+4 % ion energy [eV]
v0 = sqrt(2*K*1.602*10e-19/m);
xi = (R0 - r0) / (R0 + r0);
xi = 0.1;
Rm = Rm*1e+1;
% v = [0, sqrt(v0.^2 / Rm), sqrt(v0.^2*(1-1/Rm))]; %[x, y, z]
v = [-v0*sqrt(1/(xi+1)), v0 * sqrt(xi/(xi+1)), 1e-30]; %[x, y, z]
v_initial = v;
maxi = 1000000;
domains = linspace(1,maxi,maxi);
v_para_hist = double.empty(maxi,0);
v_perp_hist = double.empty(maxi,0);
subplot(2,2,[1,2])
trajectory1 = scatter3(double.empty(1,0), double.empty(1,0), double.empty(1,0), 'red', 'filled');
for k = 1:maxi
    xx(k) = x(1);
    yy(k) = x(2);
    zz(k) = x(3);
    vx(k) = v(1);
    vy(k) = v(2);
    vz(k) = v(3);
    [A, Bx, By, Bz] = B_field_torus(x(1), x(2), x(3), I, r0, R0, N);
    B = [Bx, By, Bz];
    v_para_hist(k) = norm(B ./ norm(B) .* dot(v, B ./ norm(B)));
    v_perp_hist(k) = norm(v - B ./ norm(B) .* dot(v, B ./ norm(B)));
    mu_hist(k) = 1/2 * m * v_perp_hist(k).^2 / norm(B);
    E_hist(k) = 1/2 * m * norm(v).^2;
    [x, v] = pusher_boris(x, v, B, m);
    if rem(k,1000) == 0 && k > 1000;
        subplot(2,2,3);
        plot(sqrt(xx.^2 + yy.^2), zz, 'r-');
        xlim([R0-r0 R0+r0])
        ylim([-r0 r0])
        subplot(2,2,4);
        hold on
        grid on
        plot(k, v_para_hist(k), 'g.', k, v_perp_hist(k), 'b.')
        xlim([0 maxi])
        ylim([-10e+1 v0+10e+1])
        legend('|v_{||}|', '|v_\perp|');
        subplot(2,2,[1,2])
        delete(trajectory1)
        trajectory1 = scatter3(xx(k-1800:k), yy(k-1800:k), zz(k-1800:k), 'red', 'filled');
        pause(0.1)
    end
    k
end

scatter3(xx, yy, zz, 'red', 'filled');
subplot(2,2,3);
plot(sqrt(xx.^2 + yy.^2), zz, 'r-');
hold on
x = linspace(0,R0+r0,100);
y = linspace(0,0,100);
z = linspace(0,0,100);
[A, Bx, By, Bz] = B_field_torus(x, y, z, I, r0, R0, N);
% plot(sqrt(x.^2 + y.^2), sqrt(Bx.^2 + By.^2)./4.1485e-04, 'b-')
plot(x, sqrt(Bx.^2 + By.^2 + Bz.^2)/0.4*r0, 'b-')

pause
figure(4)
subplot(2,2,[1,3]);
pause
for j = 1:1000:maxi
    subplot(2,2,[1,3]);
    longss = 1000;
    trajectory1 = plot(sqrt(xx(1:j).^2+yy(1:j).^2), zz(1:j), 'r.');
    axis([0 5 -1 1])
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
    hold on;
end


hold off
