clear; % Clear variables
addpath('../data')
datasetNum = 1; % CHANGE THIS VARIABLE TO CHANGE DATASET_NUM
[sampledData, sampledVicon, sampledTime, proj2Data] = init(datasetNum);
% Set initial condition
uPrev = vertcat(sampledVicon(1:9,1),zeros(6,1)); % Copy the Vicon Initial state
covarPrev = 0.01*eye(15); % Covariance constant
savedStates = zeros(15, length(sampledTime)); %Just for saving state his.
prevTime = 0;
vel = proj2Data.linearVel;
angVel2 = proj2Data.angVel;
%% Calculate Kalmann Filter
for i = 1:length(sampledTime)
    %% FILL IN THE FOR LOOP
    %ang Vel is extracted from sampledData
    angVel = sampledData(i).omg;

    %acc is extracted from sampledData
    acc = sampledData(i).acc;

    %Z_vis is calculated from 'pos' and 'pose' to use it as an input in the
    %update step
    vel_vis = transpose(vel(i,:));
    angVel2_vis = transpose(angVel2(i,:));
    Z_vis = [vel_vis;angVel2_vis];

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

plotData(savedStates, sampledTime, sampledVicon, 2, datasetNum);