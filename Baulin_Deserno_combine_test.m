R = 0.1;
d = 10*R;
phi = pi/12;
lambda = 0.01;
L0 = sin(phi)*R;
r = linspace(L0, d/2,1000);

[h,C,A] = free_shape_linear(r, L0, d, phi, lambda);

% hderiv = (C(1)./r/lambda-C(3)/lambda*besselj(1,r/lambda)-C(4)/lambda*bessely(1,r/lambda));
hderiv = (C(1)./r-C(3)/lambda*besselj(1,r/lambda)-C(4)/lambda*bessely(1,r/lambda));

% plot of shape and nanoparticle
figure();
hold on
axis equal
plot(r, h);
t = linspace(-pi/2,pi/2,100);
x = cos(t)*R;
% y = sin(t)*R+(R*cos(phi)+h(1));
y = sin(t)*R+R*cos(phi)+h(1);
plot(x,y)
% plot(x, (x-sin(phi))*tan(phi)+h(1))

% plot of derivative
% figure();
% hold on
% plot(r,hderiv)
% plot(r, tan(phi)*ones(size(r)))