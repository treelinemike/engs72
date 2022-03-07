% Cart/pendulum simulation
% Quickly updated to be driven by falling/hanging mass
% Includes some damping in pendulum joint as oscillations do tend to decay
% quickly in practice
%
% Author: M. Kokko
% Modified: 06-Mar-2022

% restart
close all; clear; clc;
if(ismac)
    % enable imagemagick and ffmpeg on mac platform 
    setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
end

% general options
anim_step = 4; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'cart_pendulum_pulled';
videoFrameRate = 30; % [frames/sec]

% simulation time parameters    
t0 = 0;         % [s] simulation start time
tf = 2;         % [s] simulation end time
dt = 0.005;      % [s] timestep size
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);

% initial conditions (state vector: [x x_dot theta theta_dot]')
x_0         = 0;  % [m]
x_dot_0     = 0;  % [m/s]
theta_0     = 0*pi/180;  % [rad]
theta_dot_0 = 0;  % [rad/s]
X0 = [x_0 x_dot_0 theta_0 theta_dot_0]'; % [m m/s rad rad/s]'
X = X0;

% system parameters
sysParams.m1 = 0.6945;    % [kg] hanging mass
sysParams.m2 = 0.4*0.4573;    % [kg] cart mass
sysParams.m3 = 0.6*0.4573;    % [kg] pendulum mass
sysParams.l  = 1;  % [m]  pendulum length
sysParams.g = 9.81;  % [m/s^2]  acceleration of gravity
sysParams.c = 0.6;   % [kg*m^2/s] damping constant for pendulum pivot

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propigate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X,opts);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% plot data
figure;

ax(1) = subplot(3,2,1);
hold on; grid on;
plot(time,data(1,:),'LineWidth',1.6);
ylabel('\bfx');

ax(2) = subplot(3,2,3);
hold on; grid on;
plot(time,data(2,:),'LineWidth',1.6);
ylabel('\bfx dot');

ax(3) = subplot(3,2,5);
hold on; grid on;
x_ddot = gradient(data(2,:),time);
plot(time,x_ddot,'LineWidth',1.6);
xlabel('\bfTime [s]');
ylabel('\bfx ddot');

ax(4) = subplot(3,2,2);
hold on; grid on;
plot(time,data(3,:),'LineWidth',1.6);
ylabel('\bftheta');

ax(5) = subplot(3,2,4);
hold on; grid on;
plot(time,data(4,:),'LineWidth',1.6);
ylabel('\bftheta dot');

ax(6) = subplot(3,2,6);
hold on; grid on;
theta_ddot = gradient(data(4,:),time);
plot(time,theta_ddot,'LineWidth',1.6);
xlabel('\bfTime [s]');
ylabel('\bftheta ddot');

linkaxes(ax,'x');

figure;
title('\bfTension in Cord');
hold on; grid on;
T = sysParams.m1*(sysParams.g - x_ddot);
plot(time,T,'LineWidth',2);
xlabel('\bfTime [s]');
ylabel('\bfTension [N]');

% % compute total energy in system
% xa = data(1,:);
% x_dot = data(2,:);
% theta = data(3,:);
% theta_dot = data(4,:);
% m1 = sysParams.m1;
% m2 = sysParams.m2;
% l = sysParams.l;
% g = sysParams.g;
% T = 0.5*m1*x_dot.^2 + 0.5*m2*(x_dot + l*theta_dot.*cos(theta)).^2 + 0.5*m2*(l*theta_dot.*sin(theta)).^2;
% V = m2*g*l*(1-cos(theta));
% E = T+V;
% 
% % plot total energy
% figure;
% hold on; grid on;
% plot(time,E-E(1),'-','LineWidth',1.6,'Color',[0 0.7 0]);
% xlabel('\bfTime [sec]');
% ylabel('\bfΔ Energy [J]');
% title('\bfEnergy Check','FontSize',12);

%% animate result in a new plot

% initialize figure
figure;
hold on; grid on;

% define cart patch geometry
cart.v = 2*[ -0.15 0;
    -0.15 0.15;
    0.15 0.15;
    0.15 0];
cart.f = [1 2 3 4 1];

% plot system
plot([-3 12],[0 0],'k-','LineWidth',2);
ph_cart  = patch('Faces',cart.f,'Vertices',nan(size(cart.v)),'FaceColor','flat','EdgeColor','k','LineWidth',2,'FaceColor',0.7*ones(1,3));
ph_pivot = plot(nan,nan,'k.','MarkerSize',40);
ph_rod   = plot(nan(1,2),nan(1,2),'k-','LineWidth',4);%,'MarkerFaceColor',[1 1 1]);
ph_bob   = plot(nan,nan,'b.','MarkerSize',40);%,'MarkerFaceColo

% finalize plot
ph_title = title('');
axis equal;
xlabel('\bfX');
ylabel('\bfY');
xlim([-1 12]);
ylim([-1.5 1]);

% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:anim_step:size(data,2)
    
    % extract state at current timestep
    y_shift    = 0.075;
    xa         = data(1,tIdx);
    theta      = data(3,tIdx);
    ya         = y_shift;
    xb         = xa+sysParams.l*sin(theta);
    yb         = ya-sysParams.l*cos(theta);
    cart.thisv = cart.v + xa*repmat([1,0],4,1);
    
    % update plot data
    ph_cart.Vertices = cart.thisv;
    ph_pivot.XData   = xa;
    ph_pivot.YData   = ya;
    ph_rod.XData     = [xa xb];
    ph_rod.YData     = [ya yb];
    ph_bob.XData     = xb;
    ph_bob.YData     = yb;
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
    system(['ffmpeg -y -r ' num2str(videoFrameRate) ' -start_number 1 -i frame%03d.png -vf "format=rgba,scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -profile:v high -pix_fmt yuv420p -g 25 -r 25 ' videoFileName '.mp4']);
    system('rm frame*.png');
end

% propagate state
function Xdot = stateProp(t,X,sysParams)

% recover parameters
m1 = sysParams.m1;
m2 = sysParams.m2;
m3 = sysParams.m3;
c = sysParams.c;

l = sysParams.l;
g = sysParams.g;

% deconstruct state vector
x         = X(1);
x_dot     = X(2);
theta     = X(3);
theta_dot = X(4);

% solve for x_ddot and theta_ddot
A = [ 0.5*m3*l*cos(theta),          (1/3)*m3*l^2;
    m1+m2+m3,   0.5*m3*l*cos(theta)];
b = [ -c*theta_dot - 0.5*m3*g*l*sin(theta);
    0.5*m3*l*(theta_dot^2)*sin(theta) + m1*g];
u = A\b;

% construct Xdot from differential equation
% note:     X    = [x x_dot theta theta_dot]'
% therefore Xdot = [x_dot x_ddot theta_dot theta_ddot]'
Xdot = zeros(4,1);
Xdot(1,:) = x_dot;
Xdot(2,:) = u(1);
Xdot(3,:) = theta_dot;
Xdot(4,:) = u(2);
end