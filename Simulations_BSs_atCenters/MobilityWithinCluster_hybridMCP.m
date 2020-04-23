function [avg_blockage_probability, avg_blockage_duration, blockage_freq] = ...
    MobilityWithinCluster_hybridMCP(BL_PPPlocs_init, BS_locs, BL_locs_initial,...
             clusterCenters,radiusCluster,Vc, Vp, L, mu, frac, simTime, tstep)
count = 0;
start_time_blockages = [0]; blockage_durations = [0];
avg_blockage_probability = 0;
avg_blockage_duration = 0;
blockage_freq = 0;

% if noBlockages(BS_locs, clusterCenters, radiusCluster, frac) == 0
%     return
% end
R=L;

intersecting_pairs = noBlockages(BS_locs, clusterCenters, radiusCluster, frac);
% size: n_BS,n_BL | tells whether a cluster circle intersects a BS link or not
dont_move_clusters = sum(intersecting_pairs,1);
% size: n_BL | tells whether how many BS links a cluster circle intersects
% if dont_move_clusters(i) = 0 then dont move cluster(i) bcoz it anyway
% never intersects any BS link

% if sum(intersecting_pairs,'all') == 0
%     return
% end

t = 0; index_step = 0;
n_BS = size(BS_locs,1);
n_BL_MCP = size(BL_locs_initial,1);
n_BL_PPP = size(BL_PPPlocs_init,1);
nsteps = ceil(simTime/tstep)+1;

blockage_timers = zeros(n_BS,(n_BL_MCP+n_BL_PPP));
BSs_blocked = zeros(n_BS,1);
BLs_timers_perBS = zeros((n_BL_MCP+n_BL_PPP),1);

%For MCP blockers
theta_MCP = -180 + rand(n_BL_MCP,1)*360;
%Vx = V.*cosd(theta_MCP); Vy = V.*sind(theta_MCP);
BL_MCP_locs_old = BL_locs_initial; BL_MCP_locs_new = BL_locs_initial;
V_vecs_MCP = zeros(n_BL_MCP,2);

%For PPP blockers
theta_PPP = -180 + rand(n_BL_PPP,1)*360;
BL_PPP_locs_old = BL_PPPlocs_init; BL_PPP_locs_new = BL_PPPlocs_init;
V_vecs_PPP = zeros(n_BL_PPP,2);

%%%%%%comment to run on HPC%%%%%%%%
%timestamps = 0:tstep:simTime;
%BL_locations = zeros(nsteps,n_BL_MCP+n_BL_PPP,2); 
%V_vectors_animation = zeros(nsteps,n_BL_MCP+n_BL_PPP,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (t < simTime)
    %%%%%%comment to run on HPC%%%%%%%%
    %BL_locations(index_step+1,:,:) =  vertcat(BL_MCP_locs_new,BL_PPP_locs_new);
    %V_vectors_animation(index_step+1,:,:) = vertcat([Vc.*cosd(theta_MCP) Vc.*sind(theta_MCP)]...
    %    ,[Vp.*cosd(theta_PPP) Vp.*sind(theta_PPP)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    V_vecs_MCP = [Vc.*cosd(theta_MCP) Vc.*sind(theta_MCP)];
    V_vecs_PPP = [Vp.*cosd(theta_PPP) Vp.*sind(theta_PPP)];
    blockage_timers = blockage_timers - tstep;
    BL_MCP_locs_new = BL_MCP_locs_old + V_vecs_MCP.*tstep;
    BL_PPP_locs_new = BL_PPP_locs_old + V_vecs_PPP.*tstep;
    for i=1:n_BL_MCP
%         if dont_move_clusters(i)==0
%             continue;
%         end
        if sum((BL_MCP_locs_new(i,:) - clusterCenters(i,:)).^2) > radiusCluster^2
            BL_MCP_locs_new(i,:) = BL_MCP_locs_old(i,:);
            theta_MCP(i) = reflectedAngle(clusterCenters(i,:), BL_MCP_locs_old(i,:),theta_MCP(i));
        end
    end
    for i=1:n_BL_PPP
        if BL_PPP_locs_new(i,1) > L || BL_PPP_locs_new(i,1) < -L
            BL_PPP_locs_new(i,1) = BL_PPP_locs_old(i,1);
            theta_PPP(i) = 180 - theta_PPP(i);
        end
        if BL_PPP_locs_new(i,2) > L || BL_PPP_locs_new(i,2) < -L
            BL_PPP_locs_new(i,2) = BL_PPP_locs_old(i,2);
            theta_PPP(i) = -theta_PPP(i);
        end
    end
    for indexBS=1:n_BS
        for indexBL_MCP=1:n_BL_MCP
            if intersecting_pairs(indexBS,indexBL_MCP) == 0
                continue;
            end
            if Intersects(BS_locs(indexBS,:)*frac,BL_MCP_locs_old(indexBL_MCP,:)...
                    ,BL_MCP_locs_new(indexBL_MCP,:))
                blockage_timers(indexBS,indexBL_MCP) = exprnd(1/mu);
            end
            
%             if blockage_timers(indexBS,indexBL) > 0
%                 blockage_instances(indexBS,index_step) = 1;
%             end
        end
        for indexBL_PPP=1:n_BL_PPP
            if Intersects(BS_locs(indexBS,:)*frac,BL_PPP_locs_old(indexBL_PPP,:)...
                    ,BL_PPP_locs_new(indexBL_PPP,:))
                blockage_timers(indexBS,n_BL_MCP+indexBL_PPP) = exprnd(1/mu);
            end
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
    
    BL_MCP_locs_old = BL_MCP_locs_new;
    BL_PPP_locs_old = BL_PPP_locs_new;
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

%RunAnimation(timestamps, clusterCenters, ...
%    radiusCluster, BS_locs, BL_locations, V_vectors_animation, 1*tstep, simTime);

end
