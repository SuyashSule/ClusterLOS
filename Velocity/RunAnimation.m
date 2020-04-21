function RunAnimation(timestamps, clusterCenters, ...
    radiusCluster, BS_locations, BL_locations, V_vectors, tstep, simTime)
hb = 1.8;
hr = 1.4;
ht = 6;
L=100;
frac = (hb-hr)/(ht-hr);
%lamT = 2e-4; %%%!!!!!!!!!!change here!!!!!!!!!!!!!!!
nB = size(BL_locations,2);
nT = size(BS_locations,1);
% rT = R*sqrt(rand(nT,1)); %location of APs
% alphaT = 2*pi*rand(nT,1);%location of APs
xT = BS_locations(:,1);%location of APs
yT = BS_locations(:,2);%location of APs
xTfrac = frac*xT; %blockage zone around UE for each APs
yTfrac = frac*yT;
v_t = 0:tstep:simTime;
for BL_index=1:nB
    vs_node(BL_index).v_x = interp1(timestamps,BL_locations(:,BL_index,1),v_t);
    vs_node(BL_index).v_y = interp1(timestamps,BL_locations(:,BL_index,2),v_t);
end

figure; hold on;
scatter(xT,yT,30,'r^')

plot([zeros(length(xT),1),xT]',[zeros(length(yT),1),yT]','g')
plot([zeros(length(xTfrac),1),xTfrac]',[zeros(length(yTfrac),1),yTfrac]',...
    'r','linewidth',1)
th=0:0.01:2*pi; xx = L*cos(th); yy = L*sin(th);
plot(xx,yy)
scatter(0,0,20,'s')

for ind=1:1:length(clusterCenters)
    %plot([zeros(length(xT),1),xT]',[zeros(length(yT),1),yT]','y')
    th=0:0.01:2*pi; 
    xx = radiusCluster*cos(th) + clusterCenters(ind,1); 
    yy = radiusCluster*sin(th) + + clusterCenters(ind,2);
    plot(xx,yy)
    scatter(0,0,20,'d')
end

for BL_index = 1:nB
    vh_node_pos(BL_index) = plot(vs_node(BL_index).v_x(1),...
        vs_node(BL_index).v_y(1),'.','color',[0.3 0.3 1],...
        'MarkerSize',10);
end
title(cat(2,'Simulation time (sec): ',num2str(simTime)));
xlabel('X (meters)');
ylabel('Y (meters)');
xlim([-100,100])
ylim([-100,100])
name=strcat('Simulation: nT=',num2str(nT),'  nB=',...
    num2str(nB), '  hT=',num2str(ht));
title(name);
ht = text(-100,90,cat(2,'Time (sec) = 0'));
%     ht = text(min(vs_node(1).v_x),max(vs_node(1).v_y),cat(2,'Time (sec) = 0'));
%     axis([min(vs_node(1).v_x) max(vs_node(1).v_x) min(vs_node(1).v_y) max(vs_node(1).v_y)]);
hold off;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!change file name
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
v = VideoWriter('simulation1.mp4');
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
v.FrameRate=10;
%     v.CompressionRatio = 3;

for timeIndex = 1:length(v_t)
    t = v_t(timeIndex);
    set(ht,'String',cat(2,'Time (sec) = ',num2str(t,4)));
    for nodeIndex = 1:nB
        set(vh_node_pos(nodeIndex),'XData',vs_node(nodeIndex).v_x(timeIndex),'YData',vs_node(nodeIndex).v_y(timeIndex));
    end
    set(gcf,'paperunits','inches')
    set(gcf,'position',[0 0 550 500])
    drawnow;
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    F(timeIndex) = getframe(ax,rect);
    
end
open(v)
writeVideo(v,F);
close(v)
end