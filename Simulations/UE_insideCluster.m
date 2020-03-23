function inside = UE_insideCluster(clusterLocations, radiusCluster)
centers = [unique(clusterLocations(:,1)) unique(clusterLocations(:,2))];
distSqFromUE = centers(:,1).^2 + centers(:,2).^2;
if sum(distSqFromUE < radiusCluster^2)~=0
    inside = 1;
else
    inside = 0;
end
end