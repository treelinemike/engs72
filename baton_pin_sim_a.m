% ENGS 72 21W Midterm Problem
%
% Simulate a baton that undergoes plastic collisions first with a fixed
% pin, and then with the ground.
%
% Author:   Mike Kokko
% Modified: 06-Mar-21

% restart
close all; clear all; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
doShowExpectedVelocities = 0;
doAnimate = 1;
anim_step = 1; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'baton_pin';
videoFrameRate = 60; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 1.94;      % [s] simulation end time
dt = 0.005;     % [s] timestep size
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);

% define parameters for cases to analyze
sysParams = [];
sysParams.m = 1;  % pendulum acceleration is mass invariant
sysParams.r1 = 0.45;
sysParams.r2 = 0.55;
sysParams.theta_pin = 30*pi/180;
sysParams.p_pin = (sysParams.r1)*[ sin(sysParams.theta_pin) cos(sysParams.theta_pin) ];

% initial mode
mode = 1;

% initial conditions (state vector: [theta theta_dot]')
theta_0     = 0.05;           % [rad]
theta_dot_0 = 0;              % [rad/s]
X0 = [theta_0 theta_dot_0]';  % [rad rad/s]'
X = X0;

% data storage
time = [t0];
data = [X0];
mode_data = [mode];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % switch modes (pivot points) and adjust for collision impulses
    switch(mode)
        case 1
            if( X(1) > sysParams.theta_pin )
                mode = 2;
                theta_dot = (sysParams.r2*(sysParams.r1+sysParams.r2)/(sysParams.r1^2+sysParams.r2^2))*X(2);
                X(2) = theta_dot;
                data(2,end) = theta_dot;
            end
        case 2
            if( (sysParams.p_pin(2) + sysParams.r2*cos(X(1))) < 0 )
                mode = 3;
                theta_dot = ((sysParams.r1)/(sysParams.r1+sysParams.r2))*X(2);
                X(2) = theta_dot;
                data(2,end) = theta_dot;
            end
    end
    
    % propigate state
    [T,X] = ode45(@(t,X) propDynamics(t,X,sysParams,mode),odeTime,X,opts);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
    mode_data(end+1) = mode;
end

%% plot trajectories

% theoretical velocities
theta_dot_2_exp = sqrt(2*9.81*(sysParams.r1+sysParams.r2)*(1-cos(sysParams.theta_pin)))/(sysParams.r1+sysParams.r2);
theta_dot_3_exp = (sysParams.r2*(sysParams.r1+sysParams.r2)/(sysParams.r1^2+sysParams.r2^2))*theta_dot_2_exp;
theta_dot_4_exp = sqrt(theta_dot_3_exp^2+2*9.81*cos(sysParams.theta_pin)*(1-(sysParams.r1/sysParams.r2))*((sysParams.r1+sysParams.r2)/(sysParams.r1^2+sysParams.r2^2)));
theta_dot_5_exp = (sysParams.r1/(sysParams.r1+sysParams.r2))*theta_dot_4_exp;
delta_E_1_2 = 0.5*sysParams.m*( ( (sysParams.r1^2+sysParams.r2^2) *theta_dot_3_exp^2) - (((sysParams.r1+sysParams.r2)^2)*theta_dot_2_exp^2) )

% initialize plot
figure;
set(gcf,'Position',[4.386000e+02 3.226000e+02 0560 2.952000e+02]);
t = tiledlayout(2,1);
t.Padding = 'none';  % 'normal', 'compact', or 'none'
t.TileSpacing = 'none';  % 'normal', 'compact', or 'none'

% angular position
ax = nexttile(1);
hold on; grid on;
plot(time,data(1,:),'-','LineWidth',2,'Color',[0 0 0.8]);
ylabel({'\bfBaton Angle [rad]',''});

% angular velocity
ax(end+1) = nexttile(2);
hold on; grid on;
plot(time,data(2,:),'-','LineWidth',2,'Color',[0 0 0.8]);
if(doShowExpectedVelocities)
    plot(time([1 end]),theta_dot_2_exp*ones(1,2),'--','Color',[0.8 0 0]);
    plot(time([1 end]),theta_dot_3_exp*ones(1,2),'--','Color',[0.8 0 0]);
    plot(time([1 end]),theta_dot_4_exp*ones(1,2),'--','Color',[0.8 0 0]);
    plot(time([1 end]),theta_dot_5_exp*ones(1,2),'--','Color',[0.8 0 0]);
end
xlabel('\bfTime [s]');
ylabel({'\bfAngular Vel. [rad/s]',''});
linkaxes(ax,'x');

%% animate result in a new plot
if(doAnimate)
    
    % initialize plot
    figure;
    hold on; grid on;
    plot([-10 10],[0 0],'-','LineWidth',6,'Color',0.7*ones(1,3));
    plot(sysParams.p_pin(1),sysParams.p_pin(2),'k.','MarkerSize',40);
    ph_rod = plot(nan(1,2),nan(1,2),'k-','LineWidth',3);
    ph_end_1 = plot(nan,nan,'.','MarkerSize',50,'Color',[ 0 0 0.8]);
    ph_end_2 = plot(nan,nan,'.','MarkerSize',50,'Color',[ 0.8 0 0]);
    axis equal;
    ph_title = title('');
    xlabel('\bfX');
    ylabel('\bfY');
    xlim([-0.5 1.0]);
    ylim([-0.25 1.25]);
   
    % animate each frame of results
    saveFrameIdx = 0;
    for tIdx = 1:anim_step:size(data,2)
        
        % extract state at current timestep
        theta = data(1,tIdx);
        theta_dot = data(2,tIdx);
        mode = mode_data(tIdx);
        
        % determine bob location
        switch(mode)
            case 1
                p_pivot = [0 0];
                p_bob_1 = p_pivot;
                p_bob_2 = (sysParams.r1+sysParams.r2)*[ sin(theta) cos(theta) ];
            case 2
                p_pivot = sysParams.p_pin;
                p_bob_1 = p_pivot - sysParams.r1*[ sin(theta) cos(theta) ];
                p_bob_2 = p_pivot + sysParams.r2*[ sin(theta) cos(theta) ];
            case 3
                p_pivot = [sysParams.p_pin(1)+sqrt(sysParams.r2^2-sysParams.p_pin(2)^2) 0];
                p_bob_1 = p_pivot + (sysParams.r1+sysParams.r2)*[-sin(theta) -cos(theta)];
                p_bob_2 = p_pivot;
        end
        
        % update graphics via blitting
        ph_rod.XData    = [p_bob_1(1) p_bob_2(1)];
        ph_rod.YData    = [p_bob_1(2) p_bob_2(2)];
        ph_end_1.XData  = p_bob_1(1);
        ph_end_1.YData  = p_bob_1(2);
        ph_end_2.XData  = p_bob_2(1);
        ph_end_2.YData  = p_bob_2(2);
        ph_title.String  = sprintf('Time: %6.3fs',time(tIdx));
        
        % finish formatting axes
        drawnow;
        
        % save frames for video if requested
        if(doMakeVideo)
            thisImgFile = sprintf('frame%03d.png',saveFrameIdx);
            saveFrameIdx = saveFrameIdx + 1;
            saveas(gcf,thisImgFile);
            system(['convert -trim ' thisImgFile ' ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
        end
    end
    
    % generate movie with ffmpeg
    if(doMakeVideo)
        system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%03d.png -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName '.mp4']);
        system('rm frame*.png');
    end
end

% propagate state
function  Xdot = propDynamics(t,X,sysParams,mode)

% deconstruct state vector
theta = X(1);
theta_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [theta      theta_dot]
% therefore Xdot = [theta_dot  theta_ddot]
Xdot = zeros(2,1);

switch(mode)
    case 1
        Xdot(1,:) = theta_dot;
        Xdot(2,:) = (9.81)*sin(theta)/(sysParams.r1+sysParams.r2);
    case 2
        Xdot(1,:) = theta_dot;
        Xdot(2,:) = (9.81)*(sysParams.r2-sysParams.r1)*sin(theta)/(sysParams.r1^2+sysParams.r2^2);
    case 3
        Xdot(1,:) = theta_dot;
        Xdot(2,:) = -(9.81)*sin(theta)/(sysParams.r1+sysParams.r2);
end
end