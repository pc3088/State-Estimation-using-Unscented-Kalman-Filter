function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state

    %'Ct' matrix defined
    Ct = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
         0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];

    %transpose(Ct) is defined
    Ct_transpose = transpose(Ct);

    %'R' matrix is defined as an identity matrix multiplied by a small
    %value
    R = eye(6)*0.01;

    %Kalman Gain is defined
    Kt = covarEst*Ct_transpose*(inv((Ct*covarEst*Ct_transpose)+R));
    
    %'uCurr' and 'covar_curr' are defined
    uCurr = uEst + Kt*(z_t - Ct*uEst);
    covar_curr = covarEst - Kt*Ct*covarEst;
    
end

