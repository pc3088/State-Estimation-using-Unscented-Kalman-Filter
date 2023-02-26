function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time
    
    %uAug is the Augmented mean calculating by concatenating uPrev and
    %zeros
    uAug = [uPrev;
            zeros(12,1)];

    %Defining Q_t which is given by Identity matrix of size 12 multiplied
    %by a small value
    Q_t = eye(12)*0.1;

    %Defining n_dash
    n_dash = 27;

    %Defining alpha
    alpha = 0.001;

    %Defining beta
    beta = 2;

    %Defining k
    k = 1;

    %Calculating lamba_dash
    lambda_dash = ((alpha^2)*(n_dash+k))-n_dash;


    root_n_dash_lambda_dash = sqrt(n_dash+lambda_dash);

    %Defining Augmented Covariance
    covarAug = [   covarPrev,  zeros(15,12);
                zeros(12,15),           Q_t];

    %Calculating square root of the Augmented Covariance by using cholesky
    %decomposition, here we are extracting the lower triangular matrix
    covarAug_sqrt = chol(covarAug,"lower");

    %Initialising sigma points as uAug in the first column
    sigma_points = uAug;
    
    %A 'for' loop is defined for Calculating the sigma points and
    %concatenating them into 'sigma_points' variable

    for i = 1:27
        X_plus  = uAug + (root_n_dash_lambda_dash*(covarAug_sqrt(:,i)));
        X_minus = uAug - (root_n_dash_lambda_dash*(covarAug_sqrt(:,i)));
        sigma_points = [sigma_points, X_plus, X_minus];
    end
    
    %Initialising a variable called 'sigma_points_output' to store the
    %output of the sigma points when passed through function 'f'
    sigma_points_output = zeros(15,55);

    %A 'for' loop is created for calculating 'sigma_points_output'
    for p = 1:55

        %Defining gravity vector 'g'
        g = [ 0; 0; -9.81];

        %Extracting roll,pitch,yaw from sigma_points
        roll  =  sigma_points(4,p);
        pitch =  sigma_points(5,p);
        yaw   =  sigma_points(6,p);

        %Defining Rotation Matrix
        R = [cos(pitch)*cos(yaw), cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw), sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch); 
             cos(pitch)*sin(yaw), cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw), cos(roll)*sin(pitch)*sin(yaw) - cos(yaw)*sin(roll);
                     -sin(pitch),                               cos(pitch)*sin(roll),                               cos(pitch)*cos(roll)];
        
        %Defining 'G' Matrix
        G = [cos(pitch)*cos(yaw), -sin(yaw), 0;
             cos(pitch)*sin(yaw),  cos(yaw), 0; 
                     -sin(pitch),         0, 1];

        %Extracting 'vm' which is the velocity
        vm =  sigma_points(7:9,p);
        
        %Extracting 'bg'
        bg = sigma_points(10:12,p);
        
        %Extracting 'ba'
        ba = sigma_points(13:15,p);
        
        %Defining elements of noise elements
        ng = sigma_points(16:18,p);
        na = sigma_points(19:21,p);
        nbg = sigma_points(22:24,p);
        nba = sigma_points(25:27,p);

        %Defing angular velocity
        wm = angVel;
        
        %Defining acceleration
        am = acc;
       
        %Defining 'sigma_points_output' variable
        sigma_points_output(:,p) = [vm;
                                    (pinv(G)*R*(wm-bg-ng));
                                    g + R*(am-ba-na);
                                    nbg;
                                    nba];

    end

    %Initialising 'Xt' for discretization 
    Xt = zeros(15,55);

    %Descritization is performed by the 'for' loop and the values are
    %stored in 'Xt'
    for i= 1:55
        Xt(:,i)= sigma_points(1:15,i) + (dt*sigma_points_output(1:15,i));
    end
    
    %W0_M,Wi_M,W0_C,Wi_C are defined for required calculations
    W0_M = (lambda_dash)/(n_dash+lambda_dash);
    Wi_M = 1/(2*(n_dash+lambda_dash));
    W0_C = ((lambda_dash)/(n_dash+lambda_dash))+(1-(alpha^2)+beta);
    Wi_C = 1/(2*(n_dash+lambda_dash));

    %mU_prev is initialized for the calculation of 'uEst'
    mU_prev = W0_M*Xt(:,1);
    
    %A 'for' loop is defined for the calculation of 'uEst'
    for i = 2:55
        mU_final = mU_prev + Wi_M*Xt(:,i);
        mU_prev  = mU_final;
    end

    uEst = mU_final;

    %A variable 'covar_previous' is initialised for the calculation of
    %'covarEst'
    covar_previous = W0_C*(Xt(:,1)-uEst)*(transpose(Xt(:,1)-uEst));
    
    %A 'for' loop is defined for the calculation of 'covarEst'
    for i = 2:55
        covar_Calculated = Wi_C*(Xt(:,i)-uEst)*(transpose(Xt(:,i)-uEst));
        covar_Final = covar_previous + covar_Calculated;
        covar_previous = covar_Final;
    end

    covarEst = covar_Final;


end