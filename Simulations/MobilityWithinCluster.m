function [avg_blockage_probability, avg_blockage_duration, blockage_freq] = ...
    MobilityWithinCluster(BS_locs, BL_locs_initial, clusterCenters, ...
    radiusCluster, V, L, mu, frac, simTime, tstep)
count = 0;
start_time_blockages = [0]; blockage_durations = [0];
avg_blockage_probability = 0;
avg_blockage_duration = 0;
blockage_freq = 0;

% if noBlockages(BS_locs, clusterCenters, radiusCluster, frac) == 0
%     return
% end

intersecting_pairs = noBlockages(BS_locs, clusterCenters, radiusCluster, frac);
dont_move_clusters = sum(intersecting_pairs,1);

if sum(intersecting_pairs,'all') == 0
    return
end

t = 0; index_step = 0;
n_BS = size(BS_locs,1);
n_BL = size(BL_locs_initial,1);
nsteps = ceil(simTime/tstep)+1;

blockage_timers = zeros(n_BS,size(BL_locs_initial,1));
BSs_blocked = zeros(n_BS,1);
BLs_timers_perBS = zeros(n_BL,1);
theta = -180 + rand(n_BL,1)*360;
Vx = V.*cosd(theta); Vy = V.*sind(theta);
BL_locs_old = BL_locs_initial; BL_locs_new = BL_locs_initial;
V_vectors = zeros(n_BL,2);

%%%%%%comment to run on HPC%%%%%%%%
% timestamps = 0:tstep:simTime;
% BL_locations = zeros(nsteps,n_BL,2); 
% V_vectors_animation = zeros(nsteps,n_BL,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (t < simTime)
    %%%%%%comment to run on HPC%%%%%%%%
%     BL_locations(index_step+1,:,:) =  BL_locs_new;
%     V_vectors_animation(index_step+1,:,:) = [V.*cosd(theta) V.*sind(theta)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_vectors = [V.*cosd(theta) V.*sind(theta)];
    blockage_timers = blockage_timers - tstep;
    BL_locs_new = BL_locs_old + V_vectors.*tstep;
    for i=1:size(BL_locs_new,1)
        if dont_move_clusters(i)==0
            continue;
        end
        if sum((BL_locs_new(i,:) - clusterCenters(i,:)).^2) > radiusCluster^2
            BL_locs_new(i,:) = BL_locs_old(i,:);
            theta(i) = reflectedAngle(clusterCenters(i,:), BL_locs_old(i,:),theta(i));
        end
    end
    for indexBS=1:n_BS
        for indexBL=1:size(BL_locs_old,1)
            if intersecting_pairs(indexBS,indexBL) == 0
                continue;
            end
            if Intersects(BS_locs(indexBS,:)*frac,BL_locs_old(indexBL,:),BL_locs_new(indexBL,:))
                blockage_timers(indexBS,indexBL) = exprnd(1/mu);
            end
            
%             if blockage_timers(indexBS,indexBL) > 0
%                 blockage_instances(indexBS,index_step) = 1;
%             end
        end
        BLs_timers_perBS = blockage_timers(indexBS,:);
        BSs_blocked(indexBS) = ~isempty(BLs_timers_perBS(BLs_timers_perBS>0));
    end    
    
    if isempty(BSs_blocked(BSs_blocked==0)) && count==0
        start_time_blockages(end) = tstep*index_step;
        start_time_blockages = [start_time_blockages 0];
        count = count + 1;
    elseif isempty(BSs_blocked(BSs_blocked==0)) && count~=0
        count = count + 1;
    elseif ~isempty(BSs_blocked(BSs_blocked==0)) && count~=0
        blockage_durations(end) = count;
        blockage_durations = [blockage_durations 0];
        count = 0;
    end
    
    BL_locs_old = BL_locs_new;
    t = t + tstep;
    index_step=index_step+1;
end

blockage_durations
start_time_blockages

sum_durs = sum(blockage_durations);

if sum_durs~=0
    avg_blockage_duration = sum_durs*tstep/(length(blockage_durations)-1);
    avg_blockage_probability = sum_durs*tstep/simTime;
    blockage_freq = (length(start_time_blockages)-1)/simTime;
end

sprintf("[avg_blockage_probability = %f, avg_blockage_duration=%f, blockage_freq=%f]"...
    ,avg_blockage_probability, avg_blockage_duration, blockage_freq)

% RunAnimation(timestamps, clusterCenters, ...
%     radiusCluster, BS_locs, BL_locations, V_vectors_animation, 10*tstep, simTime);

end
