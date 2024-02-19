% Simulate gyroscopic motion
% Simple gyroscope with length to cm = rcm and moment due to gravity
% Developed for ENGS 72, Thayer School of Engineering @ Dartmouth
% Author: Mike Kokko
% Updated: 18-Feb-2024

% Note: Damping terms can be applied to each rotation:
% - In theta: drains energy from (transient) oscilating mode
% - In phi: increases theta consistent with observed behavior
% - In psi: acts to stop spin

% restart
close all; clear; clc;

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 10;        % [s] simulation end time
dt = 0.001;     % [s] timestep size
anim_step = 10; % skip this many frames to speed up animation

% parameters of gyroscope system
params.m = 1;                  % [kg]
params.g = 9.81;               % [m/s^2]
params.I = params.m*(0.1)^2;   % [kg-m^2] common principal moment of inertia, x and y axes
params.Iz = params.m*(0.2)^2;  % [kg-m^2] 
params.rcm = 0.1;              % [m]
params.phi_damping = 0.0;      % [Nm/s]
params.theta_damping = 2.0;    % [Nm/s]   
params.psi_damping = 0.000;    % [Nm/s]

% initial conditions X0 = [phi phi_dot theta theta_dot psi psi_dot]
X0 = [0 0 pi/3 0 0 30]'; % [rad rad/s rad rad/s rad rad/s]'

% we could calcuate the phi_dot for steady precession, but it will converge
% after a short transient if we just start from an inital theta and psi_dot
% and damp out the oscillation in theta
% theta_0 = X0(3);
% psi_dot_0 = X0(6);
% phi_dot_0_candidates = roots([ (params.Iz-params.I)*cos(theta_0) params.Iz*psi_dot_0 -params.m*params.g*params.rcm]);
% X0(2) = phi_dot_0_candidates(2);
% phi_dot_0 = params.m*params.g*params.rcm/(params.Iz*psi_dot_0); % for theta = 90deg;

X = X0;

% data storage
time = [t0];
data = [X0];

%% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % net moment about pivot at this timestep
    M = zeros(3,1);
    M(1) = params.m*params.g*params.rcm*sin(X(3));

    % propagate state
    [T,X] = ode45(@(t,X) gyroStateProp(t,X,params,M),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

%% animate motion

% initialize plot
figure;
hold on; grid on; axis equal;
h_XYZ  = plotTriad(eye(4),0.02);
h_xyz  = plotTriad(nan(4),0.05);
h_body = plotTriad(nan(4),0.05);
h_title = title('');
view([120,30]);
xlim([-0.15 0.15]);
ylim([-0.15 0.15]);
zlim([-0.05 0.2]);

% step through data
for data_idx = 1:anim_step:size(data,2)

    % extract angles
    phi = data(1,data_idx);
    theta = data(3,data_idx);
    psi = data(5,data_idx);

    % assemble rotation matrices
    R1z = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
    R2x = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    R3z = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    
    % xyz frame transform
    % this is what we use as the "body frame" in gyroscopic motion
    R_xyz = R1z*R2x*eye(3);
    TF_xyz = eye(4);
    TF_xyz(1:3,1:3) = R_xyz;

    % transform for frame that fully rotates with body
    R_body = R1z*R2x*R3z*eye(3);
    TF_body = eye(4);
    TF_body(1:3,1:3) = R_body;
    TF_body(1:3,4) = R_xyz*[0; 0; params.rcm];

    % blit out updates to triads
    updateTriad(h_xyz,TF_xyz);
    updateTriad(h_body,TF_body);

    % update time and refresh display
    h_title.String = sprintf('Time: %0.2fs',time(data_idx));
    drawnow;
end


%% plot results 
figure;
ax(1) = subplot(6,1,1);
hold on; grid on;
plot(time,data(1,:));
ylabel('\bfphi [rad]');

ax(2) = subplot(6,1,2);
hold on; grid on;
plot(time,data(2,:));
ylabel('\bfphi dot [rad/s]');

ax(3) = subplot(6,1,3);
hold on; grid on;
plot(time,data(3,:));
ylabel('\bftheta [rad]');

ax(4) = subplot(6,1,4);
hold on; grid on;
plot(time,data(4,:));
ylabel('\bftheta dot [rad/s]');

ax(5) = subplot(6,1,5);
hold on; grid on;
plot(time,data(5,:));
ylabel('\bfpsi [rad]');

ax(6) = subplot(6,1,6);
hold on; grid on;
plot(time,data(6,:));
ylabel('\bfpsi dot [rad/s]');
xlabel('\bfTime [s]');

% link the x axes of each plot together
linkaxes(ax,'x');

%% state propagation function
function Xdot = gyroStateProp (t,X,params,M)

% extract parameters of interest
I = params.I;
Iz = params.Iz;

% deconstruct state vector
phi_dot = X(2);
theta = X(3);
theta_dot = X(4);
psi_dot = X(6);

% sum the moments about fixed point p (base of gyroscope
% add some damping as desired
Mx = M(1) -1*params.theta_damping*theta_dot;
My = M(2) -1*params.phi_damping*phi_dot*sin(theta);
Mz = M(3) -1*params.phi_damping*phi_dot*cos(theta) + -1*params.psi_damping*psi_dot;

% compute angular acceleration components from gyroscopic motion equations
phi_ddot =   (1/sin(theta))*(((My+Iz*theta_dot*(phi_dot*cos(theta)+psi_dot))/I)-2*phi_dot*theta_dot*cos(theta));
theta_ddot = (Mx-Iz*phi_dot*sin(theta)*(phi_dot*cos(theta)+psi_dot))/ I + (phi_dot^2)*sin(theta)*cos(theta);
psi_ddot =   (Mz/Iz)+phi_dot*theta_dot*sin(theta)-phi_ddot*cos(theta);

% construct Xdot from differential equations
% note:     X    = [phi phi_dot theta theta_dot psi psi_dot]
% therefore Xdot = [phi_dot phi_ddot theta_dot theta_ddot psi_dot psi_ddot]
Xdot = zeros(6,1);
Xdot(1,:) = phi_dot;
Xdot(2,:) = phi_ddot;
Xdot(3,:) = theta_dot;
Xdot(4,:) = theta_ddot;
Xdot(5,:) = psi_dot;
Xdot(6,:) = psi_ddot;
end