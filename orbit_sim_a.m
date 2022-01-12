% Simple Orbit Simulation
% Author: Mike Kokko
% Modified: 11-Jan-2022

% restart
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 1; % speed up animation by skipping this many frames between refreshing plot
doAnimate = 0;
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'orbit';
videoFrameRate = 400; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 200000;    % [s] simulation end time
dt = 5;         % [s] timestep size

% initial conditions

% Class Example Problem
% r_0 = 6.731e6;             % [m]
% r_dot_0 = 0;               % [m/s]
% theta_0 = 0;               % [rad]
% theta_dot_0 = 8215/(r_0);  % [rad/s]

% ISS
r_0 = 6.781e6;             % [m]
r_dot_0 = 0;               % [m/s]
theta_0 = 0;               % [rad]
theta_dot_0 = 7667/(r_0);  % [rad/s]

% PS #2
% r_0 = 2.594e7;                   % [m]
% r_dot_0 = 0;                     % [m/s]
% theta_0 = 0;                     % [rad]
% theta_dot_0 = 4166.67/(r_0);     % [rad/s]

% PS #4 Molniya Orbit
% r_0 = 7.3035e6;                  % [m]
% r_dot_0 = 0;                     % [m/s]
% theta_0 = 0;                     % [rad]
% theta_dot_0 = 7.086e10/(r_0^2);  % [rad/s]

% apply initial conditions  
X0 = [r_0 r_dot_0 theta_0 theta_dot_0]';   % [m m/s rad rad/s]'
X = X0;

% system parameters
sysParams.M = 5.976e24; % [kg] mass of Earth
sysParams.G = 6.67e-11; % [N*m^2/kg^2]  gravitational constant
sysParams.re = 6.373e6; % [m] radius of Earth
sysParams.done = false;

% data storage
time = t0;
data = X0;

% run simulation
t = t0;
while ( ~sysParams.done && (t < (tf-dt)) )

    % calculate timestep for ODE solving
    odeTime = [t t+dt];

    % check for collision with Earth or complete orbit
    if(X(1) < sysParams.re || X(3) > 2*pi)
        sysParams.done = true;
    end

    % propagate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()

    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()

    % increment to next timestep
    t = t + dt;
end

% TODO: compute total energy in system (i.e. the Hamiltonian)

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
% and other useful quantities
dvdt_cartesian = gradient(v_cartesian,time);
dvdt_polar = (r_dot.*r_ddot + r.*theta_dot.*(r_dot.*theta_dot + r.*theta_ddot))./sqrt(r_dot.^2 + (r.*theta_dot).^2);
a_radial = r_ddot;
a_centripetal = -r.*theta_dot.^2;
a_tangential = r.*theta_ddot;
a_coriolis = 2.*r_dot.*theta_dot;
a_er = a_radial + a_centripetal;
a_etheta = a_tangential+a_coriolis;

% plot results
fh = figure;
set(gcf,'Position',[0305 0233 9.568000e+02 0420]);
fh.Color = [0 0 0];
t = tiledlayout(3,5);
t.TileSpacing = 'tight';
t.Padding = 'tight';
nexttile(1,[3 2]);
hold on; grid on;  axis equal;
ah = gca;
ah.GridAlpha = 0.25;
ah.Color = 'none';%[0 0 0];
ah.XColor = [1 1 1];
ah.YColor = [1 1 1];
theta_e = 0:0.1:(2*pi+0.1);
xe = sysParams.re*cos(theta_e);
ye = sysParams.re*sin(theta_e);
% plot(xe,ye,'-','Color',[0 0 0.8],'LineWidth',1.2);
ph = patch(xe,ye,[0.4 0.4 1],'LineStyle','none','FaceAlpha',0.7);
plot(x,y,'-','Color',[1 1 0],'LineWidth',2.0);
ph_title = title('','Color',[1 1 1]);
[~,min_idx] = min(r);
[~,max_idx] = max(r);
xp = x(min_idx);
yp = y(min_idx);
xa = x(max_idx);
ya = y(max_idx);
plot([xp xa],[yp ya],':','Color',[1 1 0],'LineWidth',1.2);
plot(0,0,'.','MarkerSize',30,'Color',[0 0 0.5]);
plot(x(1),y(1),'.','MarkerSize',30,'Color',[1 0 1]);
ph_particle = plot(nan,nan,'.','MarkerSize',20,'Color',[1 1 0]);

% position plot
nexttile(3,[1,3]);
ah1 = gca;
ah1.GridAlpha = 0.25;
ah1.Color = [0 0 0];
ah1.XColor = [1 1 1];
ah1.YColor = [1 1 1];
hold on; grid on;
title('Position','Color',[1 1 1]);
ph_plist(1) = plot(time/60,r,'-','Color',[0.8 0 0],'LineWidth',2);
lh1 = legend(ph_plist,{'r'},'TextColor',[1 1 1],'EdgeColor',[1 1 1],'FontWeight','bold','Location','SouthEast');

% velocity plot
nexttile(8,[1,3]);
ah2 = gca;
ah2.GridAlpha = 0.25;
ah2.Color = [0 0 0];
ah2.XColor = [1 1 1];
ah2.YColor = [1 1 1];
hold on; grid on;
title('Speed','Color',[1 1 1]);
ph_vlist(1) = plot(time/60,v_cartesian,'-','Color',[0.8 0 0],'LineWidth',2);
ph_vlist(2) = plot(time/60,v_er,'--','Color',[0.8 .5 0],'LineWidth',2);
ph_vlist(3) = plot(time/60,v_etheta,'--','Color',[0 0.5 0.8],'LineWidth',2);
ph_v = plot(time(1)/60,nan,'.','MarkerSize',30,'Color',[0.8 0 0]);
lh2 = legend(ph_vlist,{'speed','er','etheta'},'TextColor',[1 1 1],'EdgeColor',[1 1 1],'FontWeight','bold','Location','SouthEast');

% acceleration plot
nexttile(13,[1,3]);
ah3 = gca;
ah3.GridAlpha = 0.25;
ah3.Color = [0 0 0];
ah3.XColor = [1 1 1];
ah3.YColor = [1 1 1];
hold on; grid on;
title('Acceleration','Color',[1 1 1]);
% plot(time,dvdt_polar,'--','Color',[0 0 0.8],'LineWidth',1.2);
% plot(time,a_radial,'--','Color',[0.8 .5 0],'LineWidth',2);
% plot(time,a_centripetal,'--','Color',[0 0.5 0.8],'LineWidth',2);
% plot(time,a_tangential,'--','Color',[0.5 0.8 0],'LineWidth',2);
% plot(time,a_coriolis,'--','Color',[0 0.8 0.5],'LineWidth',2);
% lh = legend('dv/dt','radial','centripetal','tangential','coriolis','TextColor',[1 1 1],'EdgeColor',[1 1 1],'FontWeight','bold');
ph_alist(1) = plot(time/60,dvdt_cartesian,'-','Color',[0.8 0 0],'LineWidth',2);
ph_alist(2) = plot(time/60,a_er,'--','Color',[0.8 .5 0],'LineWidth',2);
ph_alist(3) = plot(time/60,a_etheta,'--','Color',[0 0.5 0.8],'LineWidth',2);  % tangential and coriolis should always sum to zero in central force motion
ph_a = plot(time(1)/60,nan,'.','MarkerSize',30,'Color',[0.8 0 0]);
lh3 = legend(ph_alist,{'dv/dt','er','etheta'},'TextColor',[1 1 1],'EdgeColor',[1 1 1],'FontWeight','bold','Location','SouthEast');

% link axes
linkaxes([ah1 ah2 ah3],'x');

%% animate result in a new plot
% animate each frame of results
if(doAnimate)
    saveFrameIdx = 0;
    for tIdx = 1:anim_step:size(data,2)

        % extract phi and compute theta at this moment
        ph_particle.XData = x(tIdx);
        ph_particle.YData = y(tIdx);
        ph_v.XData = time(tIdx)/60;
        ph_v.YData = v_cartesian(tIdx);
        ph_a.XData = time(tIdx)/60;
        ph_a.YData = dvdt_cartesian(tIdx);

        % add title
        ph_title.String = sprintf('Time: %10.0f min',time(tIdx)/60);

        % draw
        drawnow;

        % save frames for video if requested
        if(doMakeVideo)
            thisImgFile = sprintf('frame%03d.png',saveFrameIdx);
            saveFrameIdx = saveFrameIdx + 1;
            exportgraphics(gcf,thisImgFile,'BackgroundColor',[0 0 0]);
            system(['convert ' thisImgFile ' -trim ' thisImgFile]);  % REQUIRES convert FROM IMAGEMAGICK!
        end
    end

    % generate movie with ffmpeg
    if(doMakeVideo)
        system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%03d.png -vf "format=rgba,scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName '.mp4']);
        system('rm frame*.png');
    end
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
end