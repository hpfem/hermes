clear all
load sol_100x100.map
[x,y,z]=deal(out(:,1),out(:,2),out(:,3));
[X,Y]=meshgrid(0:0.01:1,0:0.01:1);
Z = griddata(x,y,z,X,Y,'cubic');

figure;
surf(X,Y,Z,'FaceColor','interp','LineStyle', 'none')
lighting phong
xlabel('x')
ylabel('y')
zlabel('u(x,y)')
axis square

figure;
contourf(X,Y,Z)
xlabel('x')
ylabel('y')
axis square
