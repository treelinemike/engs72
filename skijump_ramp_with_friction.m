% Particle on a ski jump ramp (e.g. ENGS 147 ramp)
% Author: Mike Kokko
% Modified: 05-Jun-2024

% restart
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 20; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'skijump';
videoFrameRate = 100; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 2.5;       % [s] simulation end time
dt = 0.001;     % [s] timestep size

% initial conditions (state vector: [s s_dot]')
s_0 = 0;               % [m]
s_dot_0 = 0;           % [m/s]
X0 = [s_0 s_dot_0]';   % [m m/s]'
X = X0;

% system parameters
sysParams.m = 6;        % [kg] ball mass (used in energy calculation only)
sysParams.g = 9.81;     % [m/s^2]  acceleration of gravity
sysParams.mu_k = 0.1;   % coefficient of sliding friction
sysParams.mode = 0;     % 0 = initial straight ramp; 1 = curved ramp; 2 = final straight ramp
sysParams.y_min = -1.0; % where to stop the ballistic trajectory (assuming fully inelastic collisions in x and y)

% data storage
time = t0;
data_path = X0;
[x,y] = s2xy(X0(1));
data_xyz = [x;nan;y;nan];
mode = sysParams.mode;

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % set mode
    if(data_path(1,end) < 1.57 )
        sysParams.mode = 0;
        % disp('Mode 0');
    elseif(data_path(1,end) < 2.312 )
        sysParams.mode = 1;
        % disp('Mode 1');
    elseif(data_path(1,end) < 2.438)
        sysParams.mode = 2;
        % disp('Mode 2');
    elseif(sysParams.mode == 2)
        sysParams.mode = 3;
        [x,y] = s2xy(X(1));
        x_dot = X(2)*cosd(39.7);
        y_dot = X(2)*sind(39.7);
        X = [x;x_dot;y;y_dot];
        % disp('Mode 3');
    elseif(sysParams.mode == 3 && X(3) < sysParams.y_min)
        sysParams.mode = 4;
        X([2 4]) = 0;
    end

    % propagate state
    if(sysParams.mode < 3)
        % using path coordinates while on ramp
        [T,X] = ode45(@(t,X) stateProp_path(t,X,sysParams),odeTime,X);
        X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
        data_path(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
        [x,y] = s2xy(data_path(1,end));
        data_xyz(:,end+1) = [x;nan;y;nan];
    else
        % using cartesian coordinates while airborne
        [T,X] = ode45(@(t,X) stateProp_xyz(t,X,sysParams),odeTime,X);
        X = X(end, :)';
        data_xyz(:,end+1) = X;
        ds = norm(data_xyz([1 3],end)- data_xyz([1 3],end-1));
        data_path(:,end+1) = [data_path(1,end)+ds;norm(X([2,4]))];
    end

    % store results from this timestep
    time(end+1)   = T(end);
    
    mode(end+1) = sysParams.mode;
end

% TODO: compute Fn, Ff, and work done by friction

% plot results
figure;
ah(1) = subplot(4,1,1);
hold on; grid on;
plot(time,data_path(1,:),'-','LineWidth',1.6,'Color',[0 0 0.8]);
xlabel('\bfTime [sec]');
ylabel('\bfS [m]');

ah(2) = subplot(4,1,2);
hold on; grid on;
plot(time,data_path(2,:),'-','LineWidth',1.6,'Color',[0 0 0.8]);
xlabel('\bfTime [sec]');
ylabel('\bfs dot [m/s]');

ah(3) = subplot(4,1,3);
hold on; grid on;
E = zeros(size(data_path,2),1);
for data_idx = 1:size(data_xyz,2)
    s_dot = data_path(2,data_idx);
    y = data_xyz(3,data_idx);
    E(data_idx) = 9.81*y + 0.5*(s_dot)^2;
end
plot(time,E);
ylabel('\bfEnergy');
% plot(time,Fn,'-','LineWidth',1.6,'Color',[0 0 0.8]);
% plot(time,Ff,'-','LineWidth',1.6,'Color',[0.8 0 0]);
% xlabel('\bfTime [sec]');
% ylabel('\bfForce [N]');
% legend('Fn','Ff');

ah(4) = subplot(4,1,4);
hold on; grid on;
plot(time,mode);
ylabel('\bfMode');
% plot(time,U_fric,'-','LineWidth',1.6,'Color',[0.8 0 0]);
% xlabel('\bfTime [sec]');
% ylabel('\bfWork [J]');
% legend('Friction');
% 

linkaxes(ah,'x');
drawnow;

%% animate result in a new plot
fh_anim = 2;
figure(fh_anim);
% set(gcf,'Position',[1.754000e+02 2.682000e+02 0560 4.200000e+02]);
hold on; grid on;
 
%% draw ramp
s = 0:0.001:2.438;
x0 = 1.57*cosd(30)+0.6096*sind(30);
y0 = -1.57*sind(30)+0.6096*cosd(30);
x1 = x0 - 0.6096*sind(30);
y1 = y0 - 0.6096*cosd(30);
x2 = x0 + 0.6096*sind(39.7);
y2 = y0 - 0.6096*cosd(39.7);
x3 = x2 + 0.126*cosd(39.7);
y3 = y2 + 0.126*sind(39.7);
plot([0 x0 x1 x2 x3],[0 y0 y1 y2 y3],'r.','MarkerSize',20);
hold on; grid on; axis equal;
pts = [];
for s_idx = 1:length(s)
    [x,y] = s2xy(s(s_idx));
    pts(end+1,:) = [x,y];
end
plot(pts(:,1),pts(:,2),'-','LineWidth',1.6,'Color',[0 0 0]);

% draw ballistic path
air_idx = find(~isnan(data_xyz(2,:)),1,'first');
plot(data_xyz(1,air_idx:end),data_xyz(3,air_idx:end),'--','LineWidth',1.6,'Color',[0 0 0.8]);
fprintf('Launch speed: %0.2f m/s\n',data_path(2,air_idx));

% plot ball
ph_ball = plot(nan,nan,'.','MarkerSize',40,'Color',[0.8 0 0.8]);

% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:anim_step:size(data_xyz,2)
        
    % update plot
    ph_ball.XData = data_xyz(1,tIdx);
    ph_ball.YData = data_xyz(3,tIdx);

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

function [x,y] = s2xy(s)
x0 = 1.57*cosd(30)+0.6096*sind(30);
y0 = -1.57*sind(30)+0.6096*cosd(30);
x2 = x0 + 0.6096*sind(39.7);
y2 = y0 - 0.6096*cosd(39.7);
if(s < 1.57)
    x = s*cosd(30);
    y = -1*s*sind(30);
elseif(s < 2.312)
    s1 = s - 1.57;
    this_theta = (30*pi/180) - (s1/0.6096);
    x = x0 - 0.6096*sin(this_theta);
    y = y0 - 0.6096*cos(this_theta);
else
    s2 = s - 2.312;
    x = x2 + s2*cosd(39.7);
    y = y2 + s2*sind(39.7);
end
end

% propagate state
function Xdot = stateProp_path(t,X,sysParams)

% extract parameters
g = sysParams.g;
mu_k = sysParams.mu_k;
mode = sysParams.mode;

% deconstruct state vector
s     = X(1);
s_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = s_dot;
if(mode == 0)
    if(s_dot >= 0)
        Xdot(2,:) = g*(sind(30)-mu_k*cosd(30));
    else
        Xdot(2,:) = g*(sind(30)+mu_k*cosd(30));
    end
elseif(mode == 1)
    theta = (s-1.889)/0.6096;
    if(s_dot >= 0)
        Xdot(2,:) = -1*(mu_k/0.6096)*(s_dot)^2 - mu_k*g*cos(theta) - g*sin(theta);
    else
        Xdot(2,:) = (mu_k/0.6096)*(s_dot)^2 + mu_k*g*cos(theta) - g*sin(theta);
    end
elseif(mode == 2)
     if(s_dot >= 0)
        Xdot(2,:) = -g*(sind(39.7)+mu_k*cosd(39.7));
    else
        Xdot(2,:) = g*(-sind(39.7)+mu_k*cosd(39.7));
    end
else
    % TODO: add ballistic motion here
    Xdot(1,:) = 0;
    % Xdot(2,:) = 0;
end
end

% propagate state
function Xdot = stateProp_xyz(t,X,sysParams)

% extract parameters
g = sysParams.g;
mode = sysParams.mode;

% deconstruct state vector
x     = X(1);
x_dot = X(2);
y     = X(3);
y_dot = X(4);

% construct Xdot from differential equation
% note:     X    = [x x_dot y y_dot]
% therefore Xdot = [x_dot x_ddot y_dot y_ddot]
Xdot = zeros(4,1);
if(mode ~= 4)
Xdot(1,:) = x_dot;
Xdot(2,:) = 0;
Xdot(3,:) = y_dot;
Xdot(4,:) = -1*g;
end
end