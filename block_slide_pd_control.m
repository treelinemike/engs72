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
tf = 2;         % [s] simulation end time
dt = 0.015;     % [s] timestep size
opts = odeset('RelTol',1e-8,'AbsTol',1e-12);

% initial conditions (state vector: [x x_dot]')
x_0         = 0;      % [m]
x_dot_0     = 0;      % [m/s]
X0 = [x_0 x_dot_0]';  % [m m/s]'
X = X0;

% system parameters
sysParams.m   = 2;         % [kg] cart mass
sysParams.kp  = 300;       % [N/m] proportional gain
sysParams.kd  = 150;       % [Ns/m] derivative gain
sysParams.F   = 0;         % [N] computed control force to be applied
sysParams.w0  = 2*pi*1.0;  % [rad/s] driving frequency 
sysParams.xamp = 0.5;      % [m] amplitude of sinusoidal target position (assumse start from negative side)
sysParams.Fmax = 50;       % [N] saturate force at this level

% data storage
time = [t0];
data = [X0];

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % COULD SIMPLIFY THIS BY STARTING FROM ONE END AND USING 1-COS(W*T) FOR
    % POSITION AND AVOID TRANSIENT!
    x_target = sysParams.xamp*sin(sysParams.w0*odeTime(2));
    v_target = sysParams.xamp*sysParams.w0*cos(sysParams.w0*odeTime(2));  
%     x_target = sysParams.xamp*(1-cos(sysParams.w0*odeTime(2)));
%     v_target = sysParams.xamp*sysParams.w0*sin(sysParams.w0*odeTime(2));
    
    sysParams.F = sysParams.kp*(x_target-X(1)) + sysParams.kd*(v_target-X(2));
    if(abs(sysParams.F) > sysParams.Fmax)
        sysParams.F = sign(sysParams.F)*sysParams.Fmax;
    end
    
    
    % propigate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X,opts);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

%% show trajectory
x = data(1,:);
x_dot = data(2,:);

figure;
subplot(3,1,1);
hold on; grid on;
plot(time,x,'-','LineWidth',1.6);
plot(time,sysParams.xamp*sin(sysParams.w0*time),'--','Color',[0.7 0 0],'LineWidth',2);
% plot(time,sysParams.xamp*(1-cos(sysParams.w0*time)),'--','Color',[0.7 0 0],'LineWidth',2);
xlabel('\bfTime [sec]');
ylabel('\bfx [m]');

subplot(3,1,2);
hold on; grid on;
plot(time,x_dot,'-','LineWidth',1.6);
plot(time,sysParams.xamp*sysParams.w0*cos(sysParams.w0*time),'--','Color',[0.7 0 0],'LineWidth',2);
xlabel('\bfTime [sec]');
ylabel('\bfx\_dot [m]');

subplot(3,1,3);
hold on; grid on;
plot(time,gradient(x_dot,time),'-','LineWidth',1.6);
plot(time,-1*sysParams.xamp*(sysParams.w0^2)*sin(sysParams.w0*time),'--','Color',[0.7 0 0],'LineWidth',2);
xlabel('\bfTime [sec]');
ylabel('\bfx\_ddot [m]');

%% animate result in a new plot

% initialize figure
figure;
hold on; grid on;

% define cart patch geometry
cart.v = [ -0.15 0;
    -0.15 0.15;
    0.15 0.15;
    0.15 0];
cart.f = [1 2 3 4 1];

% plot system
plot([-3 3],[0 0],'k-','LineWidth',2);
ph_cart  = patch('Faces',cart.f,'Vertices',nan(size(cart.v)),'FaceColor','flat','EdgeColor','k','LineWidth',2,'FaceColor',0.7*ones(1,3));
ph_pivot = plot(nan,0.075,'k.','MarkerSize',30);

% finalize plot
ph_title = title('');
axis equal;
xlabel('\bfX');
ylabel('\bfY');
xlim([-1.0 1.0]);
ylim([-0.8 1]);
    
% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:anim_step:size(data,2)
    
    % extract state at current timestep
    
    x          = data(1,tIdx);
    cart.thisv = cart.v + x*repmat([1,0],4,1);
    
    % update plot data
    ph_cart.Vertices = cart.thisv;
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
m = sysParams.m;
F = sysParams.F;

% deconstruct state vector
x         = X(1);
x_dot     = X(2);

% construct Xdot from differential equation
% note:     X    = [x x_dot]'
% therefore Xdot = [x_dot x_ddot]'
Xdot = zeros(2,1);
Xdot(1,:) = x_dot;
Xdot(2,:) = F/m;
end