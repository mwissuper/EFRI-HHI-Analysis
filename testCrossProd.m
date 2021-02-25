% Check that taking cross product and projecting to 2D plane is equivalent
% to |r||f|sin(theta) where r and f are 2D vectors with components only in
% the plane of interest

% Care xy plane and Mz

clear; close; clc;

r = [-1 2 3];
f = [3 4 -5];

quiver(0,0,r(1),r(3)),hold on
quiver(r(1),r(3),f(1),f(3)),axis square

t1 = cross(r,f); 
ty1 = t1(2)

% theta = acos(dot(r(1:2),f(1:2)));
theta = -atan2(f(3),f(1)) + atan2(r(3),r(1));
t2 = norm(r([1 3]))*norm(f([1 3]))*sin(theta)

%% Check that y coordinate of r doesn't really matter if just care Ty

r2 = [-1 10 3];
t3 = cross(r2,f);
ty3 = t3(2)

%% 

