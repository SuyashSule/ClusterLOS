tic
aID = getenv('SLURM_ARRAY_TASK_ID')
if(isempty(aID))
  warning('aID is empty. Replacing it with 1.')  
  aID = '1'; %Runs only for first value of AP density when aID=1
end
rng('shuffle');

densityBS = [100,200,300,400,500]*10^(-6);
densityBL_PPP = [0.01,0.01,0.01];
% densityBS = 200*10^(-6);
% densityBL_PPP = 0.01;
RHO = [0.2, 0.5, 1];
L = 100;
Vp = 1; %velocity in m/s
Vc = 0.5;
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
simTime = 3*60*60; %sec Total Simulation time
% tstep = 0.0001;
tstep = 0.001; %(sec) time step
mu = 2; %Expected bloc dur =1/mu sec
R = 100; %m Radius
omegaVal = [0, pi/3]; %self blockage angle

%Parameters for Matern Cluster Process
% densityParent = [0.01/64, 0.1/400];
% radiusCluster = 1/4.0./sqrt(densityParent);
% densityDaughter = densityBL_PPP*16/pi;
numPerCl = 10;
densityParent = RHO.*densityBL_PPP/numPerCl;
radiusCluster = 10;
densityDaughter = numPerCl./pi./radiusCluster.^2;%10/pi/20/20;

%Parameter for PPP blockers
densityPPP = (1-RHO).*densityBL_PPP;

data = zeros(length(densityBS),length(densityBL_PPP),4);

for indexBS=1:length(densityBS)
     for indexBL=1:length(densityBL_PPP)
         sprintf("densityBS = %f, densityBL = %f",densityBS(indexBS),densityBL_PPP(indexBL))
         p_BS = densityBS(indexBS)/densityParent(indexBL);
         if p_BS > 1
             [numBS, BS_locs] = PPP_generate((p_BS-1)*densityBS(indexBS),L);
             p_BS = 1;
         else
             BS_locs = [];
         end
%          [numBS, BS_locs] = PPP_generate(densityBS(indexBS),L);
         [numBL_PPP, BL_PPPlocs_init] = PPP_generate(densityPPP(indexBS),L);
         [BL_locs_initial, clusterCenters, numbPointsWithinSimWindow] = ...
             MCP_generate(densityParent(indexBL),radiusCluster(indexBL),densityDaughter(indexBL),L);
         centers = unique(clusterCenters,'rows');
         BS_locs = [];
         for i=1:size(centers,1)
             if rand()<=p_BS
                 BS_locs = vertcat(BS_locs,centers(i,:));
             end
         end
         isUE_insideCluster = UE_insideCluster(clusterCenters, radiusCluster(indexBL));
         [avg_blockage_probability, avg_blockage_duration, blockage_freq] = ...
             MobilityWithinCluster_hybridMCP(BL_PPPlocs_init, BS_locs, BL_locs_initial,...
             clusterCenters,radiusCluster(indexBL),Vc, Vp, L, mu, frac, simTime, tstep);
%          outage_instances = prod(blockage_instances,1);
%          avg_blockage_probability = sum(outage_instances)/length(outage_instances);
%          blockages_durations = []; avg_blockage_duration = 0;
%          count = 0; 
%          for i=1:1:length(outage_instances)
%              if outage_instances(i)==1
%                  count = count+1;
%              elseif outage_instances(i)==0 && count~=0
%                  blockages_durations = [blockages_durations count];
%              end
%          end
%          if isempty(blockage_durations)~=0
%          avg_blockage_duration = sum(blockages_durations)/length(blockages_durations);
%          end
          data(indexBS,indexBL,:) = ...
              [isUE_insideCluster,avg_blockage_probability,avg_blockage_duration, blockage_freq];
%          data(1,1,:)
     end
end

% RunAnimation(timestamps, BS_locs, BL_locations, V_vectors, 0.1*tstep, simTime);
save(strcat('Sim_BSatCtrs_RHO_2_5_10','_',num2str(aID),'.mat'),'data')
toc
