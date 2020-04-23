function intersects = noBlockages(BS_locs, clusterCenters, radiusClusters, frac) 
%centers = [unique(clusterCenters(:,1)) unique(clusterCenters(:,2))];
BS_locs = BS_locs*frac;
intersects = zeros(size(BS_locs,1),size(clusterCenters,1));
for indexBS=1:1:size(BS_locs,1)
    for indexClusters=1:1:size(clusterCenters,1)
        e = [0,0]; l = BS_locs(indexBS,:);
        c = clusterCenters(indexClusters,:);
        r = radiusClusters;
        
        d = l - e; f = e - c;
        
        a = dot(d,d); b = 2*dot(f,d); c = dot(f,f) - r^2;
        discr = b^2 - 4*a*c;
        if discr >= 0
            discr = sqrt(discr);
            t1 = (-b - discr)/(2*a);
            t2 = (-b + discr)/(2*a);
            if t1 >= 0 && t1 <= 1
                intersects(indexBS,indexClusters) = 1;
            elseif t1<0 && t2 >0
                intersects(indexBS,indexClusters) = 1;
            end
        end
    end
end
%intersects = sum(intersects,'all');
end