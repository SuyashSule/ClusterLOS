function [avg_blockage_probability, avg_blockage_duration, blockage_freq] = ...
    MobilityWithinCluster_old(BS_locs, BL_locs_initial, clusterCenters, ...
    radiusCluster, V, L, mu, frac, simTime, tstep)
t = 0; index_step = 1;
n_BS = size(BS_locs,1);
n_BL = size(BL_locs_initial,1);
nsteps = ceil(simTime/tstep)+1;
timestamps = 0:tstep:simTime;
blockage_instances = zeros(n_BS,nsteps);
blockage_timers = zeros(n_BS,n_BL);
theta = -180 + rand(n_BL,1)*360;
Vx = V.*cosd(theta); Vy = V.*sind(theta);
BL_locs_old = BL_locs_initial; BL_locs_new = BL_locs_initial;
BL_locations = zeros(nsteps,n_BL,2);
V_vectors = zeros(nsteps,n_BL,2);
% BL_locations(1,:,:) =  BL_locs_new;
% V_vectors(1,:,:) = [Vx Vy];

while (t < simTime)
    BL_locations(index_step,:,:) =  BL_locs_new;
    V_vectors(index_step,:,:) = [V.*cosd(theta) V.*sind(theta)];
    blockage_timers = blockage_timers - tstep;
    BL_locs_new = BL_locs_old + reshape(V_vectors(index_step,:,:),n_BL,2).*tstep;
    for i=1:size(BL_locs_new,1)
        if sum((BL_locs_new(i,:) - clusterCenters(i,:)).^2) > radiusCluster^2
            BL_locs_new(i,:) = BL_locs_old(i,:);
            theta(i) = reflectedAngle(clusterCenters(i,:), BL_locs_old(i,:),theta(i));
        end
    end
    for indexBS=1:n_BS
        for indexBL=1:size(BL_locs_old,1)
            if Intersects(BS_locs(indexBS,:)*frac,BL_locs_old(indexBL,:),BL_locs_new(indexBL,:))
                blockage_timers(indexBS,indexBL) = exprnd(1/mu);
            end
            if blockage_timers(indexBS,indexBL) > 0
                blockage_instances(indexBS,index_step) = 1;
            end
        end
    end
    
    BL_locs_old = BL_locs_new;
    t = t + tstep;
    index_step=index_step+1;
end

outage_instances = prod(blockage_instances,1);
avg_blockage_probability = sum(outage_instances)/length(outage_instances);
blockage_durations = [0]; avg_blockage_duration = 0;

count = 0;
for i=1:1:length(outage_instances)
    if outage_instances(i)==1
        count = count+1;
    elseif outage_instances(i)==0 && count~=0
        blockage_durations = [blockage_durations count];
        count = 0;
    end
end
blockage_durations
%if isempty(blockage_durations)~=0
    avg_blockage_duration = sum(blockage_durations)*tstep/length(blockage_durations);
%end

blockage_freq = 0;
for i=1:1:length(outage_instances)-1
    if outage_instances(i+1) - outage_instances(i)==1
        blockage_freq = blockage_freq+1;
    end
end
blockage_freq = blockage_freq/simTime;

sprintf("[avg_blockage_probability = %f, avg_blockage_duration=%f, blockage_freq=%f]"...
    ,avg_blockage_probability, avg_blockage_duration, blockage_freq)

RunAnimation(timestamps, clusterCenters, ...
    radiusCluster, BS_locs, BL_locations, V_vectors, 10*tstep, simTime);

end
