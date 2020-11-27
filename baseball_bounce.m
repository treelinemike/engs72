% Baseball Bounce Simulation
% For introducing basics of using ode45 for solving simple ODEs numerically

% restart
close all; clear; clc;

% general options
animate_step = 10; % speed up animation by skipping this many frames between refreshing plot
doMakeVideo = 1; % set to 1 to produce a video file; requires imagemagick ('convert') and ffmpeg
videoFileName = 'baseball_drop';
videoFrameRate = 100; % [frames/sec]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 1.15;      % [s] simulation end time
dt = 0.001;     % [s] timestep size

% initial conditions (state vector: [y y_dot]')
y_0     = 1;         % [m]
y_dot_0 = 0;         % [m/s]
X0 = [y_0 y_dot_0]'; % [m m/s]'
X = X0;

% system parameters
sysParams.m = 1;    % [kg]     problem is mass independent!
sysParams.g = 9.81; % [m/s^2]  acceleration of gravity
sysParams.e = 0.45; % coefficient of restitution

% data storage
time = t0;
data = X0;

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % apply bounce if needed
    if( X(1) < 0 && X(2) < 0 )
        X(2) = -1*sysParams.e*X(2);
    end
    
    % propagate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

% compute total energy in system
y = data(1,:);
y_dot = data(2,:);
E = 0.5*sysParams.m*y_dot.^2 + sysParams.m*sysParams.g*y;

% plot results
figure;
ah(1) = subplot(3,1,1);
hold on; grid on;
plot(time,data(1,:),'b-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfPosition [m]');

ah(2) = subplot(3,1,2);
hold on; grid on;
plot(time,data(2,:),'r-','LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfVelocity [m/s]');

ah(3) = subplot(3,1,3);
hold on; grid on;
plot(time,E,'-','LineWidth',1.6,'Color',[0 0.7 0]);
xlabel('\bfTime [sec]');
ylabel('\bfEnergy [J]');

linkaxes(ah,'x');

%% Animate result in a new plot
figure;
hold on; grid on;
xlim_vals = 0.5*0.618*[-1 1];

% draw ground
patch(xlim_vals([1 1 2 2 ]),[0 -0.1 -0.1 0],[0.6 0.6 0.6],'EdgeColor','none');

% draw ball
ph_ball = plot(0,NaN,'o','MarkerSize',12,'LineWidth',4,'MarkerFaceColor',[1 1 1],'Color',[0.8 0 0]);

% format plot
axis equal;
ylabel('\bfHeight [m]');
xlim(xlim_vals);
ylim([-0.1 1]);
ax = gca;
ax.XAxis.Visible = 'off';

% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:animate_step:size(data,2)
    
    % extract state at current timestep
    y = data(1,tIdx);
    
    % update ball position 
    ph_ball.YData = y;
    
    % finish formatting axes
    title(sprintf('Time: %6.3fs',time(tIdx)));
	drawnow;
%     pause(0.001);

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

% extract parameters
g = sysParams.g;

% deconstruct state vector
y     = X(1);
y_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = y_dot;
Xdot(2,:) = -1*g;  

end