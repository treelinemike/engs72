% Orbit Simulation
% Author: Mike Kokko
% Modified: 27-Dec-2021

% restart
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 3; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'orbit';
videoFrameRate = 100; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 30000;     % [s] simulation end time
dt = 1;         % [s] timestep size

% initial conditions
r_0 = 6.731e6;             % [m]
r_dot_0 = 1000;               % [m/s]
theta_0 = 0;               % [rad]
theta_dot_0 = 8215.8/r_0;  % [rad/s]
X0 = [r_0 r_dot_0 theta_0 theta_dot_0]';   % [m m/s rad rad/s]'
X = X0;

% system parameters
sysParams.M = 5.976e24; % [kg] mass of Earth
sysParams.G = 6.67e-11; % [N*m^2/kg^2]  gravitational constant
sysParams.done = false;
sysParams.re = 6.373e6; % [m] radius of Earth

% data storage
time = t0;
data = X0;

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    if(X(1) < sysParams.re)
       sysParams.done = true;
    end
    % propagate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% TODO: compute total energy in system

% extract data and compute derivatives
r = data(1,:);
r_dot = data(2,:);
r_ddot = gradient(r_dot,time);
theta = data(3,:);
theta_dot = data(4,:);
theta_ddot = gradient(theta_dot,time);

x = r.*cos(theta);
y = r.*sin(theta);
dxdt = gradient(x,time);
dydt = gradient(y,time);
v_cartesian = sqrt(dxdt.^2 + dydt.^2);
v_er = r_dot;
v_etheta = r.*theta_dot;
v_polar = sqrt(v_er.^2 + v_etheta.^2);  % MAGNITUDE of velocity vector (i.e. speed)

% COMPUTE TIME RATE OF CHANGE OF SPEED (not magnitude of acceleration)
dvdt_cartesian = gradient(v_cartesian,time);
% a_er = r_ddot;
% a_etheta = r.*theta_ddot + 2.*r_dot.*theta_dot;
% dvdt_polar = sqrt(a_er.^2 + a_etheta.^2);
dvdt_polar = (r_dot.*r_ddot + r.*theta_dot.*(r_dot.*theta_dot + r.*theta_ddot))./sqrt(r_dot.^2 + (r.*theta_dot).^2);

% plot results
figure;
hold on; grid on; axis equal;
plot(0,0,'.','MarkerSize',20,'Color','k');
theta_e = 0:0.1:(2*pi+0.1);
xe = sysParams.re*cos(theta_e);
ye = sysParams.re*sin(theta_e);
plot(xe,ye,'-','Color',[0 0 0.8],'LineWidth',1.2);
plot(x,y,'-','Color',[0.8 0 0.8],'LineWidth',2.0);
[~,min_idx] = min(r);
[~,max_idx] = max(r);
xp = x(min_idx);
yp = y(min_idx);
xa = x(max_idx);
ya = y(max_idx);
plot([xp xa],[yp ya],'--','Color',[0.8 0 0.8],'LineWidth',1.2);

figure;
subplot(2,1,1);
hold on; grid on;
plot(time,v_cartesian,'-','Color',[0.8 0 0],'LineWidth',1.2);
plot(time,v_polar,'--','Color',[0 0 0.8],'LineWidth',1.2);


subplot(2,1,2);
hold on; grid on;
plot(time,dvdt_cartesian,'-','Color',[0.8 0 0],'LineWidth',1.2);
plot(time,dvdt_polar,'--','Color',[0 0 0.8],'LineWidth',1.2);


%% animate result in a new plot
fh_anim = 2;
figure(fh_anim);
% set(gcf,'Position',[1.754000e+02 2.682000e+02 0560 4.200000e+02]);
hold on; grid on;
 
% draw bowl
theta_bowl = pi:0.01:2*pi;
x_bowl = sysParams.R*cos(theta_bowl);
y_bowl = sysParams.R*sin(theta_bowl);
plot(x_bowl,y_bowl,'-','LineWidth',6,'Color',[0 0 0]);
xlim_vals = sysParams.R*[-1 1];
ylim_vals = sysParams.R*[-1.2 0.5];
xlim(xlim_vals);
ylim(ylim_vals);
axis equal;

% template for ball
ball_segs = 4;
theta_ball = 0:0.01:2*pi;
x_ball_template = sysParams.r*cos(theta_ball);
y_ball_template = sysParams.r*sin(theta_ball);
N_pts_ball_seg = floor(length(theta_ball)/ball_segs);
ball_patch_v_template = [0 0; x_ball_template(1:N_pts_ball_seg)' y_ball_template(1:N_pts_ball_seg)'; 0 0];
ball_patch_v = [];
ball_patch_f = [];

% extend template around ball for patch 'stripes'
for ballSegIdx = 1:2:ball_segs
    theta_rotate = (ballSegIdx-1)*(2*pi/ball_segs);
    R = [cos(theta_rotate) -sin(theta_rotate);sin(theta_rotate) cos(theta_rotate)];
    startIdx = size(ball_patch_v,1)+1;
    endIdx = startIdx + size(ball_patch_v_template,1) -1;
    ball_patch_v(startIdx:endIdx,:) = ball_patch_v_template*R;
    ball_patch_f = [ball_patch_f; (startIdx + (1:size(ball_patch_v_template,1)) -1 )];
    
end
ball_patch_v_template = ball_patch_v;

% apply template ball at intial position
this_phi = X0(1);
this_theta = (-1*(sysParams.R-sysParams.r)/sysParams.r)*this_phi;
r_ctr = (sysParams.R-sysParams.r)*[sin(this_phi) -cos(this_phi)];
R = [cos(this_theta) sin(this_theta); -sin(this_theta) cos(this_theta)];
x_ball = x_ball_template + r_ctr(1);
y_ball = y_ball_template + r_ctr(2);
ball_patch_v = ball_patch_v_template*R+r_ctr;
ph_ball_patch = patch('Vertices',ball_patch_v,'Faces',ball_patch_f,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
ph_ball_od = plot(x_ball,y_ball,'-','LineWidth',4','Color',[0.8 0 0]);


% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:anim_step:size(data,2)
    
    % extract phi and compute theta at this moment
    this_phi = data(1,tIdx);
    this_theta = (-1*(sysParams.R-sysParams.r)/sysParams.r)*this_phi;
    
    % apply ball template at this moment
    r_ctr = (sysParams.R-sysParams.r)*[sin(this_phi) -cos(this_phi)];
    R = [cos(this_theta) sin(this_theta); -sin(this_theta) cos(this_theta)];
    x_ball = x_ball_template + r_ctr(1);
    y_ball = y_ball_template + r_ctr(2);
    ball_patch_v = ball_patch_v_template*R+r_ctr;
    ph_ball_patch.Vertices = ball_patch_v;
    
    % update plot
    ph_ball_od.XData = x_ball;
    ph_ball_od.YData = y_ball;
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
G = sysParams.G;
M = sysParams.M;

% deconstruct state vector
r         = X(1);
r_dot     = X(2);
theta_dot = X(4);

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = zeros(4,1);
Xdot(1,:) = r_dot;
Xdot(2,:) = r*theta_dot^2 - (G*M)/(r^2);
Xdot(3,:) = theta_dot;
Xdot(4,:) = -1*(2*r_dot*theta_dot)/ r;
Xdot = Xdot*(~sysParams.done); % squash this if we're done
end