% pushing a block simulation

% restart
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 1; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'block_slide_pd_control';
videoFrameRate = 30; % [frames/sec]

% simulation time parameters    
t0 = 0;         % [s] simulation start time
tf = 4;         % [s] simulation end time
dt = 0.015;     % [s] timestep size
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);

% initial conditions (state vector: [x x_dot theta theta_dot]')
theta_0         = 0;  % [rad]
theta_dot_0     = 0;  % [rad/s]
X0 = [theta_0 theta_dot_0]'; % [rad rad/s]'
X = X0;

% system parameters
sysParams.md     = 2.0;  % [kg] mass of block
sysParams.rd     = 0.05;  % [m] half height of block
sysParams.l      = 0.125; % [m] half length of block
sysParams.Izz_cm = (1/12)*sysParams.md*((2*sysParams.rd)^2+(2*sysParams.l)^2); % [kg*m^2] block moment of inertia about its CM
sysParams.g      = 9.81;  % [m/s^2] acceleration of gravity
sysParams.k      = 500;   % [N/m] spring constant
sysParams.c      = 0.2;     % [Ns/m] damping constant
sysParams.w0     = 2*pi*0.5; % [rad/s] forcing frequency
sysParams.xamp   = 0.5;   % [m] amplitude of displacement vibration
sysParams.x      = 0;     % storage for current value of linear acceleration
sysParams.x_ddot = 0;     % storage for current value of linear acceleration

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];

    % compute current position and acceleration of sled (and thus pivot point)
    sysParams.x = sysParams.xamp*(1-cos(sysParams.w0*odeTime(2)));
    sysParams.x_ddot = sysParams.xamp*(sysParams.w0^2)*cos(sysParams.w0*odeTime(2));  
    
    % propigate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X,opts);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

%% compute total energy in system
theta = data(1,:);
x_dot = data(2,:);
% theta = data(3,:);
% theta_dot = data(4,:);
% ma = sysParams.ma;
% mb = sysParams.mb;
% l = sysParams.l;
% g = sysParams.g;
% T = 0.5*ma*x_dot.^2 + 0.5*mb*(x_dot + l*theta_dot.*cos(theta)).^2 + 0.5*mb*(l*theta_dot.*sin(theta)).^2;
% V = mb*g*l*(1-cos(theta));
% E = T+V;
% 
% % % plot trajectory
figure;
subplot(4,1,1);
hold on; grid on;
plot(time,data(1,:)*180/pi,'-','Color',[0.7 0 0],'LineWidth',2);
xlabel('\bfTime [sec]');
ylabel('\bfBlock Angle [deg]');

% plot(time,sysParams.xamp*(1-cos(sysParams.w0*time)),'--','Color',[0.7 0 0],'LineWidth',2);
% 
% subplot(4,1,2);
% hold on; grid on;
% plot(time,x_dot,'-','LineWidth',1.6);
% plot(time,sysParams.xamp*sysParams.w0*cos(sysParams.w0*time),'--','Color',[0.7 0 0],'LineWidth',2);
% 
% subplot(4,1,3);
% hold on; grid on;
% plot(time,gradient(x_dot,time),'-','LineWidth',1.6);
% plot(time,-1*sysParams.xamp*(sysParams.w0^2)*sin(sysParams.w0*time),'--','Color',[0.7 0 0],'LineWidth',2);
% 
% 

% plot(time,E-E(1),'-','LineWidth',1.6,'Color',[0 0.7 0]);
% xlabel('\bfTime [sec]');
% ylabel('\bfΔ Energy [J]');
% title('\bfEnergy Check','FontSize',12);
% 
%% animate result in a new plot

% initialize figure
figure;
set(gcf,'Position',[1.602000e+02 3.026000e+02 1.002400e+03 0420]);
hold on; grid on;

% define cart patch geometry
block.v = [ -sysParams.l 0;
    -sysParams.l 2*sysParams.rd;
    sysParams.l 2*sysParams.rd;
    sysParams.l 0];
block.f = [1 2 3 4 1];

% plot system
plot([-3 3],[0 0],'k-','LineWidth',2);
ph_block  = patch('Faces',block.f,'Vertices',nan(size(block.v)),'FaceColor','flat','EdgeColor','k','LineWidth',2,'FaceColor',0.7*ones(1,3));
ph_pivot = plot(nan,0,'k.','MarkerSize',30);

% finalize plot
ph_title = title('');
axis equal;
xlabel('\bfX');
ylabel('\bfY');
xlim([-0.5 1.5]);
ylim([-0.2 .2]);
    
% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:anim_step:size(data,2)
    
    % extract state at current timestep
    t = time(tIdx);
    
    theta          = data(1,tIdx);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    x = sysParams.xamp*(1-cos(sysParams.w0*t));
    block.thisv = block.v*R + x*repmat([1,0],4,1);
    
    % update plot data
    ph_block.Vertices = block.thisv;
    ph_pivot.XData   = x;
    ph_title.String  = sprintf('Time: %6.3fs',time(tIdx));
    
    % redraw plot
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

% propagate state
function Xdot = stateProp(t,X,sysParams)

% recover parameters
md      = sysParams.md;
x_ddot  = sysParams.x_ddot;
Izz_cm  = sysParams.Izz_cm;
g       = sysParams.g;
k       = sysParams.k;
c       = sysParams.c;
rd      = sysParams.rd;
l       = sysParams.l;

% deconstruct state vector
theta     = X(1);
theta_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [x x_dot theta theta_dot]'
% therefore Xdot = [x_dot x_ddot theta_dot theta_ddot]'
Xdot = zeros(2,1);
Xdot(1,:) = theta_dot;
Xdot(2,:) = (md*rd*g*sin(theta)-2*(k*l^2 + md*rd^2*theta_dot^2)*sin(theta)*cos(theta) -md*rd*x_ddot*cos(theta) - 2*c*theta_dot*(l*cos(theta))^2) /(Izz_cm - md*rd^2);
end