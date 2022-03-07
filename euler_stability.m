% Simulate Euler's equations of motion for a rectangular board
% Developed for ENGS 72, Thayer School of Engineering @ Dartmouth
% Author: Mike Kokko
% Updated: 18-Feb-2021

% restart
close all; clear; clc;

% simulation time parameters
t0 = 0;          % [s] simulation start time
tf = 20;         % [s] simulation end time
dt = 0.001;      % [s] timestep size
anim_step = 50;  % skip this many frames to speed up animation

% physical parameters of board
rho = 0.4e3;                         % [kg/m^3] density
l_x = 0.25;                          % [m] thickness along x
l_y = 0.10;                          % [m] thickness along y
l_z = 0.01;                          % [m] thickness along z
V   = l_x*l_y*l_z;                   % [m^3] volume
m   = rho*V;                         % [kg] mass
Ixx = (1/12)*m*(l_y^2+l_z^2);        % moment of inertia about body x axis
Iyy = (1/12)*m*(l_x^2+l_z^2);        % moment of inertia about body y axis
Izz = (1/12)*m*(l_x^2+l_y^2);        % moment of inertia about body z axis
Icm  = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];  % inertia matrix taken at CM about principal axes

% initial conditions X0 = [theta_x_0 theta_y_0 theta_z_0 omega_x_0 omega_y_0 omega_z_0]
X0 = [0 0 0 0 2 0.1]'; % [rad rad rad rad/s rad/s rad/s]'
X = X0;

% data storage
time = [t0];
data = [X0];

%% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    % propagate state
    [T,X] = ode45(@(t,X) simpleEulerSimStateProp(t,X,Icm),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

%% develop and display animation of motion
% define Cartesian frames
h = 0.9*max([l_x, l_y, l_z]); % scaling parameter
triad_XYZ = h*[0 0 0; 1 0 0; 0 0 0; 0 1 0; 0 0 0; 0 0 1];
triad_xyz_template = .8*triad_XYZ;

% assemble model of board in body coordintes
board.v = [ ...
    l_x/2, -l_y/2,  l_z/2;
    l_x/2,  l_y/2,  l_z/2;
    l_x/2,  l_y/2, -l_z/2;
    l_x/2, -l_y/2, -l_z/2;
    -l_x/2, -l_y/2,  l_z/2;
    -l_x/2,  l_y/2,  l_z/2;
    -l_x/2,  l_y/2, -l_z/2;
    -l_x/2, -l_y/2, -l_z/2 ];
board.f = [ ...
    1, 2, 3, 4, 1;
    4, 3, 7, 8, 4;
    8, 7, 6, 5, 8;
    5, 1, 2, 6, 5;
    5, 1, 4, 8, 5;
    6, 2, 3, 7, 6 ];
board.c = repmat([0.2 0.4 0.6 0.8 0 0]',1,3); % face colors

% keep track of body axes to transform H back
R = eye(3);

% initialize figure
figure;
hold on; grid on;

% plot each structure and grab handles for blitting
ph.XYZ = plot3(triad_XYZ(:,1),triad_XYZ(:,2),triad_XYZ(:,3),'LineWidth',4,'Color','k');
ph.xyz_x = plot3(nan(1,2),nan(1,2),nan(1,2),'LineWidth',3,'Color',[0.7 0 0]);
ph.xyz_y = plot3(nan(1,2),nan(1,2),nan(1,2),'LineWidth',3,'Color',[0 0.7 0]);
ph.xyz_z = plot3(nan(1,2),nan(1,2),nan(1,2),'LineWidth',3,'Color',[0 0 0.7]);

% plot angular velocity and angular momentum for blitting
ph.omega = plot3(nan(1,2),nan(1,2),nan(1,2),':','LineWidth',3','Color',[1 0 1]);
ph.H = plot3(nan(1,2),nan(1,2),nan(1,2),':','LineWidth',3','Color',[0 1 1]);

% plot patch body for blitting
ph.patch = patch('Faces',board.f,'Vertices',nan(size(board.v)),'FaceColor','flat','EdgeColor','none','LineWidth',1,'FaceVertexCData',board.c);

% add additional plot features
th = title('');
legend([ph.omega,ph.H],{'Angular Velocity','Angular Momentum'},'Location','southoutside','AutoUpdate','off');
xlabel('\bfx');
ylabel('\bfy');
zlabel('\bfz');
view([145,30]);

% animate each frame of the results
for tIdx = 2:size(data,2)
    
    % compute rotation between body coordinates and inertial space
    d_theta_x = data(1,tIdx) - data(1,tIdx-1);
    d_theta_y = data(2,tIdx) - data(2,tIdx-1);
    d_theta_z = data(3,tIdx) - data(3,tIdx-1);
    xRot = [1 0 0 ; 0 cos(d_theta_x) -sin(d_theta_x); 0 sin(d_theta_x) cos(d_theta_x)];
    yRot = [cos(d_theta_y) 0 sin(d_theta_y); 0 1 0; -sin(d_theta_y) 0 cos(d_theta_y)];
    zRot = [cos(d_theta_z) -sin(d_theta_z) 0 ; sin(d_theta_z) cos(d_theta_z) 0; 0 0 1];
    R  = R*(xRot*yRot*zRot*eye(3));     % incremental rotations so order shouldn't matter
    
    % skip frames if desired to speed up animation
    % don't do this in the for loop b/c need to update rotation at each step
    if( mod(tIdx-2,anim_step) == 0 )
        % transform point cloud to correct location in inertial space
        board.vrot = board.v*R';
        triad_xyz = triad_xyz_template*R';
        
        % compute angular velocity in terms of the inertial basis (XYZ)
        omega = data(4:6,tIdx);
        omega_XYZ = R*omega;
        
        % compute angular momentum vector in terms of inertial basis (XYZ)
        % THIS SHOULD STAY CONSTANT (no external moments)
        Hcm_xyz = Icm*omega;
        Hcm_XYZ = R*Hcm_xyz;
        
        % blit the body xyz frame
        ph.xyz_x.XData = [triad_xyz(1,1) triad_xyz(2,1)];
        ph.xyz_x.YData = [triad_xyz(1,2) triad_xyz(2,2)];
        ph.xyz_x.ZData = [triad_xyz(1,3) triad_xyz(2,3)];
        ph.xyz_y.XData = [triad_xyz(3,1) triad_xyz(4,1)];
        ph.xyz_y.YData = [triad_xyz(3,2) triad_xyz(4,2)];
        ph.xyz_y.ZData = [triad_xyz(3,3) triad_xyz(4,3)];
        ph.xyz_z.XData = [triad_xyz(5,1) triad_xyz(6,1)];
        ph.xyz_z.YData = [triad_xyz(5,2) triad_xyz(6,2)];
        ph.xyz_z.ZData = [triad_xyz(5,3) triad_xyz(6,3)];
        
        % normalize and plot angular velocity and momentum
        omega_norm = 2.6*omega_XYZ/norm(omega_XYZ);
        Hcm_norm = 2.6*Hcm_XYZ/norm(Hcm_XYZ);
        ph.omega.XData = [0 omega_norm(1)];
        ph.omega.YData = [0 omega_norm(2)];
        ph.omega.ZData = [0 omega_norm(3)];
        ph.H.XData = [0 Hcm_norm(1)];
        ph.H.YData = [0 Hcm_norm(2)];
        ph.H.ZData = [0 Hcm_norm(3)];
                
        % plot board as patch object
        ph.patch.Vertices = board.vrot;
        
        % finish formatting axes
        axis equal;
        xlim([-h h]);
        ylim([-h h]);
        zlim([-h h]);       
        th.String = sprintf('Euler Motion Sim (%6.3fs)',time(tIdx));
        drawnow;
    end
end

%% plot angular velocity trajectories
figure;
ax = subplot(3,1,1);
hold on; grid on;
plot(time,data(4,:),'-','Color',[0.8 0 0],'LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfw_x [rad/s]');

ax(end+1) = subplot(3,1,2);
hold on; grid on;
plot(time,data(5,:),'-','Color',[0 0.8 0],'LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfw_y [rad/s]');

ax(end+1) = subplot(3,1,3);
hold on; grid on;
plot(time,data(6,:),'-','Color',[0 0 0.8],'LineWidth',1.6);
xlabel('\bfTime [sec]');
ylabel('\bfw_z [rad/s]');

% set axis limits
linkaxes(ax,'xy');
xlim(time([1 end]));
ylim([ min([0, min(min(data(4:6,:)))]) max([0, max(max(data(4:6,:)))]) ]);

%% state propagation function
function Xdot = simpleEulerSimStateProp(t,X,Icm)

% recover moments of inertia
Ixx = Icm(1,1);
Iyy = Icm(2,2);
Izz = Icm(3,3);

% deconstruct state vector
theta_x = X(1);
theta_y = X(2);
theta_z = X(3);
omega_x = X(4);
omega_y = X(5);
omega_z = X(6);

% construct Xdot from differential equation
% note:     X    = [theta_x theta_y theta_z omega_x omega_y omega_z]
% therefore Xdot = [omega_x omega_y omega_z alpha_x alpha_y alpha_z]
Xdot = zeros(6,1);
Xdot(1,:) = omega_x;
Xdot(2,:) = omega_y;
Xdot(3,:) = omega_z;
Xdot(4,:) = (omega_y*omega_z*(Iyy-Izz))/Ixx;
Xdot(5,:) = (omega_z*omega_x*(Izz-Ixx))/Iyy;
Xdot(6,:) = (omega_x*omega_y*(Ixx-Iyy))/Izz;
end