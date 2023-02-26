clear; % Clear variables
addpath('../data')
datasetNum = 1; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime,proj2Data] = init(datasetNum);
Z = sampledVicon(1:6,:);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.1*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0; %last time step in real time
pos = proj2Data.position;
pose = proj2Data.angle;

for i = 1:length(sampledTime)
    %% Fill in the FOR LOOP
    %ang Vel is extracted from sampledData
    angVel = sampledData(i).omg;

    %acc is extracted from sampledData
    acc = sampledData(i).acc;

    %Z_vis is calculated from 'pos' and 'pose' to use it as an input in the
    %update step
    pos_vis = transpose(pos(i,:));
    pose_vis = transpose(pose(i,:));
    Z_vis = [pos_vis;pose_vis];

    %dt is defined
    dt = sampledTime(i)- prevTime;

    %Prediction step
    [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt);

    %Update Step
    [uCurr,covar_curr] = upd_step(Z_vis,covarEst,uEst);

    %The calculated 'uCurr' from the update step is saved into
    %'savedStates'
    savedStates(:,i) = uCurr;
    
    %'uPrev' and 'covarPrev' are re-defined as 'uCurr' and 'covar_curr' for the next loop
    uPrev = uCurr;
    covarPrev = covar_curr;

    %prevTime is initialised again
    prevTime = sampledData(i).t;


end

plotData(savedStates, sampledTime, sampledVicon, 1, datasetNum);
