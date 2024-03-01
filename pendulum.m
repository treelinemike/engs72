% Simulate and plot phase portraits for the simple 1DOF pendulum
% Demo for ENGS 72: Dynamics, Thayer School of Engineering @ Dartmouth College
% Author: Mike Kokko
% Modified: 01-Mar-2024

% restart
close all; clear all; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 10; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'pendulum';
videoFrameRate = 20; % [frames/sec]
doShowVectorField = true;
doShowLagrangian = false;

% simulation time parameters
t0 = 0;      % [s] simulation start time
tf = 10;     % [s] simulation end time
dt = 0.01;   % [s] timestep size
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

% define parameters for cases to analyze
sysParams = [];
sysParams.m = 1;  % pendulum acceleration is mass invariant
sysParams.l = 6;
sysParams.c = 0.0; % try 0 or 0.6
sysParams.g = 9.81; % [m/s^2]
sysParams.T_boost = 0;

% initial conditions (state vector: [theta theta_dot]')
theta_0     = -120*pi/180;    % [rad]
theta_dot_0 = 0;              % [rad/s]
X0 = [theta_0 theta_dot_0]';  % [rad rad/s]'
X = X0;

% data storage
time = [t0];
data = [X0];

boost_timer = 0;
boost_dir = 1;

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
%     if( (abs(X(1)) < (pi/6)) )
%         sysParams.T_boost = sign(X(2))*1.2;
%     end

    % propigate state
    [T,X] = ode45(@(t,X) propDynamics(t,X,sysParams),odeTime,X,opts);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% % check energy
% theta = data(1,:);
% theta_dot = data(2,:);
% E = sysParams.m*sysParams.g*sysParams.l*(1-cos(theta)) + 0.5*sysParams.m*(sysParams.l*theta_dot).^2;
% figure;
% hold on; grid on;
% plot(time,E-E(1));
% xlabel('\bfTime');
% ylabel('\bfΔ Energy [J]');

% configure phase portrait grid
[x1,x2] = meshgrid(-6.5:0.5:6.5,-4:0.5:4);

% compute phase portrait
x1dot = x2;
x2dot = -1*(sysParams.g/sysParams.l)*sin(x1)-sysParams.c*x2;
% x2dot = -1*(sysParams.g/sysParams.l)*sin(x1)-sysParams.c*x2 + 1*(pi-x1) + (sysParams.g/sysParams.l)*sin(x1) - 2.5*x2;
% x2dot = -1*(sysParams.g/sysParams.l)*sin(x1)-sysParams.c*x2 + 2*(sysParams.g/sysParams.l)*sin(x1);

% plot phase portrait
figure;
figpos = get(gcf,'Position');
deltafig = figpos(3:4) - [834.4000 558.4000];
set(gcf,'Position',[figpos(1:2)+deltafig figpos(3:4)-deltafig] );
t = tiledlayout(3,3);
t.Padding = 'none';  % 'normal', 'compact', or 'none'
t.TileSpacing = 'none';  % 'normal', 'compact', or 'none'


% plot phase plane
nexttile(2,[2 2]);
hold on; grid on;
if(doShowVectorField)
    quiver(x1,x2,x1dot,x2dot,1.5,'LineWidth',2,'Color',0.6*ones(1,3));
end
xlabel('\bfAngular Position [rad]');
ylabel('\bfAngular Velocity [rad/s]');
title('Pendulum Phase Portrait');
axis equal;

% plot fixed points in phase plane
fixedPtLocs = [-2*pi 0; -pi 0; 0 0; pi 0; 2*pi 0]';
plot(fixedPtLocs(1,1:2:size(fixedPtLocs,2)),fixedPtLocs(2,1:2:size(fixedPtLocs,2)),'ko','MarkerSize',15,'LineWidth',4,'MarkerFace','k');
plot(fixedPtLocs(1,2:2:size(fixedPtLocs,2)),fixedPtLocs(2,2:2:size(fixedPtLocs,2)),'ko','MarkerSize',15,'LineWidth',4);

% plot lagrangian
T = 0.5*(sysParams.m*sysParams.l)^2*x2.^2;
V = sysParams.m*sysParams.g*sysParams.l*(1-cos(x1));
L = T-V;
if(doShowLagrangian)
    contour(x1,x2,L);
end

%% animate result

% plot phase space trajectory
x1traj = data(1,:)';
x2traj = data(2,:)';
plot(X0(1),X0(2),'.','Color',[0 0.8 0],'MarkerSize',40);
ph_phase_traj = plot(nan(1,length(time)),nan(1,length(time)),'-','Color',[0 0.8 0],'LineWidth',4);
set(gca,'XLim',[-6 6]);
set(gca,'YLim',[-4 4]);

% plot time series
nexttile(7,[1 3]);
hold on; grid on;
plot(time,0*ones(1,length(time)),'k-');
ph_tseries_pos = plot(time,nan(1,length(time)),'-','Color',[0 0 0.8],'LineWidth',2);
ph_tseries_vel = plot(time,nan(1,length(time)),'-','Color',[0.8 0 0],'LineWidth',2);
legend([ph_tseries_pos,ph_tseries_vel],{'Position [rad]','Velocity [rad/s]'},'Location','Northeast');
xlabel('\bfTime [s]');
set(gca,'XLim',[t0 tf]);
set(gca,'YLim',[min([x1traj;x2traj]),max([x1traj;x2traj])]);

% plot pendulum
nexttile(1,[2,1]);
hold on; grid on;
plot(0,0,'.','MarkerSize',40,'Color',0.6*ones(1,3));
ph_bob_line = plot(nan(1,2),nan(1,2),'-','LineWidth',6,'Color',0.6*ones(1,3));
ph_bob_end = plot(nan,nan,'o','MarkerSize',20,'LineWidth',6,'MarkerFaceColor',0.3*ones(1,3),'Color',0.3*ones(1,3));
axis equal;
xlim(1.2*sysParams.l*[-1 1]);
ylim(1.6*sysParams.l*[-1 1]);
axh = gca;
axh.XAxis.Visible = 'off';
axh.YAxis.Visible = 'off';
ph_time_text = text(0,8,sprintf('Time: %6.3fs',0),'HorizontalAlignment','center','FontWeight','bold','FontSize',14);

% animate each frame of results
saveFrameIdx = 0;
thisDataIdx = 1;
for tIdx = 1:anim_step:size(data,2)
    
    % extract state at current timestep
    theta = data(1,tIdx);
    theta_dot = data(2,tIdx);
   
    % recover length
    l = sysParams.l;
    
    % determine bob location
    bob_x =  l*sin(theta);
    bob_y = -l*cos(theta);
    
    % update figure via blitting
    ph_phase_traj.XData(thisDataIdx)  = x1traj(tIdx);
    ph_phase_traj.YData(thisDataIdx)  = x2traj(tIdx);
    ph_tseries_pos.XData(thisDataIdx) = time(tIdx);
    ph_tseries_pos.YData(thisDataIdx) = x1traj(tIdx);
    ph_tseries_vel.XData(thisDataIdx) = time(tIdx);
    ph_tseries_vel.YData(thisDataIdx) = x2traj(tIdx);
    ph_bob_line.XData    = [0 bob_x];
    ph_bob_line.YData    = [0 bob_y];
    ph_bob_end.XData     = bob_x;
    ph_bob_end.YData     = bob_y;
    ph_time_text.String  = sprintf('Time: %6.3fs',time(tIdx));
	thisDataIdx = thisDataIdx + 1;
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
function  Xdot = propDynamics(t,X,sysParams)

% recover paramters
m = sysParams.m;  % note: pendulum acceleration is mass invariant
l = sysParams.l;
c = sysParams.c;
g = sysParams.g;
T_boost = sysParams.T_boost;

% deconstruct state vector
theta = X(1);
theta_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [theta      theta_dot]
% therefore Xdot = [theta_dot  theta_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = theta_dot;
Xdot(2,:) = -1*(g/l)*sin(theta)-c*theta_dot + T_boost;
% Xdot(2,:) = -1*(g/l)*sin(theta)-c*theta_dot +  (g/l)*sin(theta) + 1*(pi-theta) - 2.5*theta_dot;


end
