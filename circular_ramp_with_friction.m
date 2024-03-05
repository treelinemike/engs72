% Particle on circular path with friction
% Author: Mike Kokko
% Modified: 04-Mar-2024

% restart
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 50; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'circular_ramp';
videoFrameRate = 100; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 0.3;       % [s] simulation end time
dt = 0.0001;    % [s] timestep size

% initial conditions (state vector: [s s_dot]')
s_0 = 0.4189;            % [m]
s_dot_0 = 2.0;           % [m/s]
X0 = [s_0 s_dot_0]';     % [m m/s]'
X = X0;

% system parameters
sysParams.m = 6;        % [kg] ball mass (used in energy calculation only)
sysParams.g = 9.81;     % [m/s^2]  acceleration of gravity
sysParams.R = 1.2;      % [m] radius of bowl, measured via CMM 09-Mar-21 (Kori Jevsevar)
sysParams.mu_k = 0.3;   % coefficient of sliding friction
sysParams.mode = 0;     % 0 = on ramp; 1 = ballistic; 2 = no motion

% data storage
time = t0;
data = X0;

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % check Fn
    Fn = sysParams.m*sysParams.g*cos((1/sysParams.R)*X(1)) - (sysParams.m/sysParams.R)*(X(2))^2;
    %%if((sysParams.mode == 0) && abs(X(1)-0.6981*sysParams.R) < 0.001)  % to check work done by friction
    if((sysParams.mode == 0) && (Fn < 0))
        sysParams.mode = 2;
        theta_end_deg = (X(1)/sysParams.R)*180/pi;
        fprintf('Done at theta = %04.2f°\n',theta_end_deg);
    end

    % propagate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% compute Fn, Ff, and work done by friction
s = data(1,:);
s_dot = data(2,:);
theta = (1/sysParams.R)*s;
Fn = sysParams.m*sysParams.g*cos(theta)-(sysParams.m/sysParams.R)*(s_dot.^2);
Ff = sysParams.mu_k*Fn;
U_fric = cumtrapz(s,Ff);

% plot results
figure;
ah(1) = subplot(4,1,1);
hold on; grid on;
plot(time,data(1,:),'-','LineWidth',1.6,'Color',[0 0 0.8]);
xlabel('\bfTime [sec]');
ylabel('\bfS [m]');

ah(2) = subplot(4,1,2);
hold on; grid on;
plot(time,data(2,:),'-','LineWidth',1.6,'Color',[0 0 0.8]);
xlabel('\bfTime [sec]');
ylabel('\bfs dot [m/s]');

ah(3) = subplot(4,1,3);
hold on; grid on;
plot(time,Fn,'-','LineWidth',1.6,'Color',[0 0 0.8]);
plot(time,Ff,'-','LineWidth',1.6,'Color',[0.8 0 0]);
xlabel('\bfTime [sec]');
ylabel('\bfForce [N]');
legend('Fn','Ff');

ah(4) = subplot(4,1,4);
hold on; grid on;
plot(time,U_fric,'-','LineWidth',1.6,'Color',[0.8 0 0]);
xlabel('\bfTime [sec]');
ylabel('\bfWork [J]');
legend('Friction');


linkaxes(ah,'x');
drawnow;

%% animate result in a new plot
fh_anim = 2;
figure(fh_anim);
% set(gcf,'Position',[1.754000e+02 2.682000e+02 0560 4.200000e+02]);
hold on; grid on;
 
% draw ramp
theta_bowl = 0:0.01:(pi/2);
x_bowl = sysParams.R*sin(theta_bowl);
y_bowl = sysParams.R*cos(theta_bowl);
plot(x_bowl,y_bowl,'-','LineWidth',6,'Color',[0 0 0]);
xlim_vals = sysParams.R*[-0.1 1.25];
ylim_vals = sysParams.R*[-0.1 1.25];
xlim(xlim_vals);
ylim(ylim_vals);
axis equal;

% plot ball
ph_ball = plot(nan,nan,'.','MarkerSize',60,'Color',[0.8 0 0]);

% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:anim_step:size(data,2)
        
    % compute bal position
    this_theta = data(1,tIdx)/sysParams.R;
    x_ball = sysParams.R*sin(this_theta);
    y_ball = sysParams.R*cos(this_theta);
    
    % update plot
    ph_ball.XData = x_ball;
    ph_ball.YData = y_ball;

    figure(fh_anim); % this is silly, needed to work around an issue with plot losing focus in MATLAB Online
    title(sprintf('Time: %6.3fs',time(tIdx)));
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
    system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%03d.png -vf "format=rgba,scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName '.mp4']);
    system('rm frame*.png');
end


% propagate state
function Xdot = stateProp(t,X,sysParams)

% extract parameters
g = sysParams.g;
R = sysParams.R;
mu_k = sysParams.mu_k;
mode = sysParams.mode;

% deconstruct state vector
s     = X(1);
s_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = zeros(2,1);
if(mode == 0)
    Xdot(1,:) = s_dot;
    Xdot(2,:) = (mu_k/R)*(s_dot)^2 + g*sin((1/R)*s) - mu_k*g*cos((1/R)*s);
else
    Xdot(1,:) = 0;
    Xdot(2,:) = 0;
end
end