function [locations, clusterLocations, numbPointsWithinSimWindow] = MCP_generate(lambdaParent,radiusCluster,lambdaD,L)
%Simulate a Matern cluster point process on a rectangle.
%Author: H. Paul Keeler, 2018.

%Simulation window parameters
xMin=-L;
xMax=L;
yMin=-L;
yMax=L;

%Parameters for the parent and daughter point processes
%lambdaPPP=0.1;
%lambdaParent=lambdaPPP/400;%density of parent Poisson point process
%radiusCluster=1/4.0/sqrt(lambdaParent);%radius of cluster disk (for daughter points)
%lambdaD = lambdaPPP*16/pi; %density of daughter Poisson point process
meanPointsDaughters=pi*radiusCluster^2*lambdaD;%mean number of points in each cluster

%Extended simulation windows parameters
rExt=radiusCluster; %extension parameter -- use cluster radius
xMinExt=xMin-rExt;
xMaxExt=xMax+rExt;
yMinExt=yMin-rExt;
yMaxExt=yMax+rExt;
%rectangle dimensions
xDeltaExt=xMaxExt-xMinExt;
yDeltaExt=yMaxExt-yMinExt;
areaTotalExt=xDeltaExt*yDeltaExt; %area of extended rectangle

%Simulate Poisson point process for the parents
numbPointsParent=poissrnd(areaTotalExt*lambdaParent,1,1);%Poisson number 
%x and y coordinates of Poisson points for the parent
xxParent=xMinExt+xDeltaExt*rand(numbPointsParent,1);
yyParent=yMinExt+yDeltaExt*rand(numbPointsParent,1);

%Simulate Poisson point process for the daughters (ie final poiint process)
numbPointsDaughter=poissrnd(meanPointsDaughters,numbPointsParent,1); 
numbPoints=sum(numbPointsDaughter); %total number of points

%Generate the (relative) locations in polar coordinates by 
%simulating independent variables.
theta=2*pi*rand(numbPoints,1); %angular coordinates 
rho=radiusCluster*sqrt(rand(numbPoints,1)); %radial coordinates 

%Convert from polar to Cartesian coordinates
[xx0,yy0]=pol2cart(theta,rho);

%replicate parent points (ie centres of disks/clusters) 
xxP=repelem(xxParent,numbPointsDaughter);
yyP=repelem(yyParent,numbPointsDaughter);
%translate points (ie parents points are the centres of cluster disks)
xx=xxP(:)+xx0;
yy=yyP(:)+yy0;

%thin points if outside the simulation window
%booleInside=((xx>=xMin)&(xx<=xMax)&(yy>=yMin)&(yy<=yMax));
booleInside=(xx.^2 + yy.^2 <= L^2);
%retain points inside simulation window
xx_retained=xx(booleInside); 
yy_retained=yy(booleInside); 
xxP_retained = xxP(booleInside);
yyP_retained = yyP(booleInside);

numbPointsWithinSimWindow = size(xx_retained,1);

locations = [xx_retained yy_retained];
clusterLocations = [xxP_retained yyP_retained];

%Plotting
%figure;
%scatter(xx,yy,'.');
%shg;
end