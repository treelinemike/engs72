% Ball in Bowl Simulation
% Using nonlinear differential equation
% Activity for ENGS 72: Dynamics, Thayer School of Engineering at Dartmouth
% Author: Mike Kokko
% Modified: 01-Mar-2021

% restart
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 3; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 1; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'ball_in_bowl';
videoFrameRate = 100; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 1.0;       % [s] simulation end time
dt = 0.001;     % [s] timestep size

% initial conditions (state vector: [phi phi_dot]')
phi_0     = 45*pi/180;   % [rad]
phi_dot_0 = 0;           % [rad/s]
X0 = [phi_0 phi_dot_0]'; % [rad rad/s]'
X = X0;

% system parameters
sysParams.m = 0.003384; % [kg] ball mass
sysParams.g = 9.81;     % [m/s^2]  acceleration of gravity
sysParams.R = 0.05808;  % [m] radius of bowl, measured via CMM 09-Mar-21 (Kori Jevsevar)
sysParams.r = 0.0095;   % [m] radius of ball

% data storage
time = t0;
data = X0;

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
        
    % propagate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% compute total energy in system
phi = data(1,:);
phi_dot = data(2,:);
theta_dot = (-1*(sysParams.R-sysParams.r)/sysParams.r)*phi_dot;
Izz_cm = (2/5)*sysParams.m*sysParams.r^2;
E = 0.5*sysParams.m*((sysParams.R-sysParams.r)*phi_dot).^2 + 0.5*Izz_cm*(theta_dot.^2) + sysParams.m*sysParams.g*(sysParams.R-sysParams.r)*(1-cos(phi));

% plot results
figure;
ah(1) = subplot(3,1,1);
hold on; grid on;
plot(time,data(1,:),'b-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfPosition [rad]');

ah(2) = subplot(3,1,2);
hold on; grid on;
plot(time,data(2,:),'r-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfSpeed [rad/s]');

ah(3) = subplot(3,1,3);
hold on; grid on;
plot(time,E-E(1),'-','LineWidth',1.6,'Color',[0 0.7 0]);
xlabel('\bfTime [sec]');
ylabel('\bfΔ Energy [J]');

linkaxes(ah,'x');
drawnow;

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

% % save data
% save('ball_05.mat','data','time');
% 
% % restart
% close all; clear; clc;
% 
% % format figure
% figure;
% set(gcf,'Position',[3.162000e+02 3.274000e+02 7.752000e+02 3.416000e+02]);
% hold on; grid on;
% title('\bfBall in Bowl Natural Period');
% xlabel('\bfTime [s]');
% ylabel('\bfAngle [rad]');
% 
% % load and plot data
% load('ball_90.mat');
% ph1 = plot(time,data(1,:),'-','LineWidth',1.6,'Color',[0.8 0 0]);
% tIdx = 1+find(data(1,2:end) > data(1,2),1,'first');
% datatip(ph1,time(tIdx),data(1,tIdx),'Location','southeast');
% 
% load('ball_45.mat');
% ph2 = plot(time,data(1,:),'-','LineWidth',1.6,'Color',[0 0.8 0]);
% tIdx = 1+find(data(1,2:end) > data(1,2),1,'first');
% datatip(ph2,time(tIdx),data(1,tIdx),'Location','southeast');
% 
% load('ball_05.mat');
% ph3 = plot(time,data(1,:),'-','LineWidth',1.6,'Color',[0 0 0.8]);
% tIdx = 1+find(data(1,2:end) > data(1,2),1,'first');
% datatip(ph3,time(tIdx),data(1,tIdx),'Location','southeast');
% 
% % add legend
% legend('90°','45°','5°');


% propagate state
function Xdot = stateProp(t,X,sysParams)

% extract parameters
g = sysParams.g;
R = sysParams.R;
r = sysParams.r;

% deconstruct state vector
phi     = X(1);
phi_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = phi_dot;
Xdot(2,:) = -1*g*sin(phi)/((7/5)*(R-r));  

end