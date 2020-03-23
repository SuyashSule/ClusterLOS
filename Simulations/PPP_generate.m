function [numbPoints, locations] = PPP_generate(density,L)
R = L;
area = pi*R^2;
lambda = density; 
numbPoints=poissrnd(area*lambda);
while numbPoints==0
    numbPoints=poissrnd(area*lambda);
end
r = R*sqrt(rand(numbPoints,1));
theta = 2*pi*(rand(numbPoints,1)-0.5);
x = r.*cos(theta);
y = r.*sin(theta);

locations = [x y];
 
%Plotting
% figure;
% scatter(xx,yy,'.');
% xlabel('x');ylabel('y');
end