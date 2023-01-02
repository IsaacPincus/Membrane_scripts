dx1 = 1;
dx2 = 1;

figure();
hold on
axis equal
line([-dx1/2 dx1/2 dx1/2 -dx1/2 -dx1/2], [-dx2/2 -dx2/2 dx2/2 dx2/2 -dx2/2],...
    'lineStyle', '--', 'color', 'k')

dy1 = 2;
dy2 = 0.5;

line([-dy1/2 dy1/2 dy1/2 -dy1/2 -dy1/2], [-dy2/2 -dy2/2 dy2/2 dy2/2 -dy2/2],...
    'lineStyle', '-', 'color', 'k')

xlim([-max(dx1,dy1), max(dx1,dy1)])
ylim([-max(dx2,dy2), max(dx2,dy2)])

lam1 = dy1/dx1;
lam2 = dy2/dx2;