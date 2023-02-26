function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    
    %Defining n
    n = 15;

    %Defining alpha
    alpha = 0.001;

    %Defining beta
    beta = 2;

    %Defining k
    k = 1;

    %Calculating lamba
    lambda = ((alpha^2)*(n+k))-n;
    
    root_n_lambda = sqrt(n+lambda);

    %W0_M,Wi_M,W0_C,Wi_C are defined for required calculations
    W0_M = (lambda)/(n+lambda);
    Wi_M = 1/(2*(n+lambda));
    W0_C = ((lambda)/(n+lambda))+(1-(alpha^2)+beta);
    Wi_C = 1/(2*(n+lambda));

    %Defining the additive noise and multiplying by a small value
    Rt = eye(3)*0.01;

    %Defining the required translation and rotation matrices between camera
    %and body frames, these were same values calculated in Project 2
    BTC = [-0.04; 0; -0.03];

    BRC = [ cosd(45), -sind(45), 0;
           -sind(45), -cosd(45), 0;
                   0,         0, -1];

    CRB = transpose(BRC);

    %Calculating Homogeneous Matrices
    BHC = [BRC, BTC;
            0,0,0,1];
    CHB = inv(BHC);

    %Initialising sigma points in the measurement model
    sigma_points_m = uEst;
    
    %Calculating Cholesky decomposition of covariance
    sqrt_covariance_m = chol(covarEst, "lower");

    %'for' loop for computing sigma points
    for i = 1:n

        sigma_points_m_plus = uEst + root_n_lambda*sqrt_covariance_m(:,i);
        sigma_points_m_minus = uEst - root_n_lambda*sqrt_covariance_m(:,i);
        sigma_points_m = [sigma_points_m, sigma_points_m_plus, sigma_points_m_minus];
    end

    %'For' loop for calculating 'Zt'
    for i = 1:(2*n+1)
        roll = sigma_points_m(4,i);
        pitch = sigma_points_m(5,i);
        yaw = sigma_points_m(6,i);

        Rot_matrix = [cos(pitch)*cos(yaw),cos(yaw)*sin(pitch)*sin(roll)-cos(roll)*sin(yaw),sin(roll)*sin(yaw)+cos(roll)*cos(yaw)*sin(pitch);
                      cos(pitch)*sin(yaw),cos(roll)*cos(yaw)+sin(pitch)*sin(roll)*sin(yaw),cos(roll)*sin(pitch)*sin(yaw)-cos(yaw)*sin(roll);
                      -sin(pitch),cos(pitch)*sin(roll),cos(pitch)*cos(roll)];
        
        %Formula for calculating Zt
        Zt(:,i) = (BRC*transpose(Rot_matrix)*(sigma_points_m(7:9,i))) - (BRC*skew(CHB(1:3,4))*CRB*z_t(4:6));
    end

    %'For' loop for calculating 'z'
    for i = 1 : (2*n + 1)
        if i == 1
          z = W0_M*Zt(:,i);
        else
        z = z + (Wi_M * Zt(:,i));
        end
    end

    %'For' loop for calculating Ct and St
    for i = 1 : (2*n + 1)
        if i == 1
          Ct = W0_C*(sigma_points_m(1:15 ,i) - uEst) * transpose(Zt(:, i) - z);
          St = W0_C * (Zt(:, i) - z) * transpose(Zt(:, i) - z);
        else
          Ct = Ct + (Wi_C * (sigma_points_m(1:15 ,i) - uEst) * transpose(Zt(:, i) - z));
          St = St + (Wi_C * (Zt(:, i) - z) * transpose(Zt(:, i) - z));
        end
    end

    %Calculating 'St' by adding 'Qt'
    St = St + Rt;

    %calculating kalman gain
    K_t = Ct/(St);  

    %Calculating uCurr and covar_curr
    uCurr = uEst + (K_t * (z_t(1:3) - z));
    covar_curr = covarEst - (K_t * St * transpose(K_t));
    



end

