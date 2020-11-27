% Flying cable simulation
% Feynman exercise #10.15

% restart
close all; clear; clc;

% general options
anim_step = 20;  % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 0; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'flying_cable';
videoFrameRate = 50; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 1.5;       % [s] simulation end time
dt = 0.001;     % [s] timestep size

% initial conditions (state vector: [y y_dot]')
y_0     = 0.005;          % [m]       slight perturbation!
y_dot_0 = 0;              % [m/s]
X0      = [y_0 y_dot_0]'; % [m m/s]'
X       = X0;             % [m m/s]'

% system parameters
sysParams.L = 1;          % [m]       total length of cable
sysParams.g = 9.81;       % [m/s^2]   acceleration of gravity
sysParams.rho = 1;        % [kg/m]    rope mass per unit length

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

% compute total energy in system (i.e. the Hamiltonian)
% depends upon whether the cable is unwrapping around pulley (E1)
% or falling freely (E2)
y = data(1,:);
y_dot = data(2,:);
m1 = (sysParams.L/2)-y;
m2 = (sysParams.L/2)+y;
m = m1+m2; % (== sysParams.L !)
E1 = 0.5*(sysParams.rho*sysParams.L)*y_dot.^2 + sysParams.rho*sysParams.g*(0.25*sysParams.L^2-y.^2);
E2 = 0.5*(sysParams.rho*sysParams.L)*y_dot.^2 + sysParams.rho*sysParams.g*(0.5*sysParams.L-y);
Emask = (y <= 0.5*sysParams.L);
E = [E1(Emask)';E2(~Emask)'];

% determine position and velocity at time when cable leaves pulley
targetIdx = find(y >= sysParams.L/2,1,'first');
fprintf('Cable leaves pulley at t = %5.3fs with v = %5.3fm/s\n',time(targetIdx),y_dot(targetIdx));

% plot results
fid_plots = 1;
figure(fid_plots);
ah(1) = subplot(3,1,1);
hold on; grid on;
plot(time([1,end]),y(targetIdx)*ones(1,2),'-','LineWidth',1,'Color',[0.8 0.0 0.8]);
plot(time,data(1,:),'-','LineWidth',1.6,'Color',[0.0 0.0 0.8]);
xlabel('\bfTime [sec]');
ylabel('\bfDisplacement [m]');

ah(2) = subplot(3,1,2);
hold on; grid on;
plot(time([1,end]),y_dot(targetIdx)*ones(1,2),'-','LineWidth',1,'Color',[0.8 0.0 0.8]);
plot(time,data(2,:),'r-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfVelocity [m/s]');

ah(3) = subplot(3,1,3);
hold on; grid on;
plot(time,E,'-','LineWidth',1.6,'Color',[0 0.7 0]);
xlabel('\bfTime [sec]');
ylabel('\bfEnergy [J]');

linkaxes(ah,'x');
drawnow;

%% Animate result in a new plot
figure;
hold on; grid on;

% parameters of the pulley and cable wrapped around it (wrap length not
% accounted for in simulation!)
r_pulley = 0.1;
theta = 0:0.01:pi;
x_circ = r_pulley*cos(theta);
y_circ = r_pulley*sin(theta);

% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:anim_step:size(data,2)
    
    % extract state at current timestep
    y = data(1,tIdx);
    y_dot = data(2,tIdx);
   
    % clear axes and start plotting the current frame
    cla;
    plot(0,0,'.','MarkerSize',50,'Color',[0 0 0]);
    if(y <= sysParams.L/2)
        plot(-r_pulley*ones(1,2), [y-sysParams.L/2 0],'-','LineWidth',3,'Color',[0.8 0 0]);
        plot(x_circ,y_circ,'-','LineWidth',3,'Color',[0.8 0 0]);
        plot(r_pulley*ones(1,2), [0 -(sysParams.L/2+y)],'-','LineWidth',3,'Color',[0.8 0 0]);
    else
        plot(r_pulley*ones(1,2),(-y+(sysParams.L/2))+[0 -sysParams.L],'-','LineWidth',3,'Color',[0.8 0 0]);
    end
    
    % finish formatting axes
    axis equal;
    xlabel('\bfX');
    ylabel('\bfY');
    xlim([-0.4 0.4]);
    ylim([-2.5 2*r_pulley]);
    title(sprintf('Time: %6.3fs',time(tIdx)));
 	drawnow;
%     pause(0.1);
    
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
L = sysParams.L;
g = sysParams.g;

% deconstruct state vector
y     = X(1);
y_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = y_dot;

% equation of motion changes based on whether
% cable is unwinding or falling freely
if(y <= L/2)
    Xdot(2,:) = 2*y*g/L;
else
    Xdot(2,:) =  g;
end
end