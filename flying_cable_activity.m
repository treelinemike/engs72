% Flying cable simulation activity for ENGS 72
% Feynman exercise #10.15

% restart
close all; clear; clc;

% general options
anim_step = 20;  % speed up animation by skipping this many frames between refreshing plot

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 1.5;       % [s] simulation end time
dt = 0.001;     % [s] timestep size

% initial conditions (state vector: [s s_dot]')
s_0     = 0.01;          % [m]       slight perturbation!
s_dot_0 = 0;              % [m/s]
X0      = [s_0 s_dot_0]'; % [m m/s]'
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
s = data(1,:);
s_dot = data(2,:);
m1 = (sysParams.L/2)-s;
m2 = (sysParams.L/2)+s;
m = m1+m2; % (== sysParams.L !)
E1 = 0.5*(sysParams.rho*sysParams.L)*s_dot.^2 + sysParams.rho*sysParams.g*(0.25*sysParams.L^2-s.^2);
E2 = 0.5*(sysParams.rho*sysParams.L)*s_dot.^2 + sysParams.rho*sysParams.g*(0.5*sysParams.L-s);
Emask = (s <= 0.5*sysParams.L);
E = [E1(Emask)';E2(~Emask)'];

% determine position and velocity at time when cable leaves pulley
targetIdx = find(s >= sysParams.L/2,1,'first');
fprintf('Cable leaves pulley at t = %5.3fs with v = %5.3fm/s\n',time(targetIdx),s_dot(targetIdx));

% plot results
fid_plots = 1;
figure(fid_plots);
ah(1) = subplot(3,1,1);
hold on; grid on;
plot(time([1,end]),s(targetIdx)*ones(1,2),'-','LineWidth',1,'Color',[0.8 0.0 0.8]);
plot(time,data(1,:),'-','LineWidth',1.6,'Color',[0.0 0.0 0.8]);
t_end = find( time <= sqrt(2*data(1,end)/sysParams.g),1,'last'); 
plot(time(1:t_end),0.5*9.81*time(1:t_end).^2,'--','LineWidth',1.6,'Color',[0.0 0.0 0]);
xlabel('\bfTime [sec]');
ylabel('\bfDisplacement [m]');

ah(2) = subplot(3,1,2);
hold on; grid on;
plot(time([1,end]),s_dot(targetIdx)*ones(1,2),'-','LineWidth',1,'Color',[0.8 0.0 0.8]);
plot(time,data(2,:),'r-','LineWidth',1.6);
t_end = find( time <= (data(2,end)/sysParams.g),1,'last'); 
plot(time(1:t_end),9.81*time(1:t_end),'--','LineWidth',1.6,'Color',[0.0 0.0 0]);
xlabel('\bfTime [sec]');
ylabel('\bfSpeed [m/s]');

ah(3) = subplot(3,1,3);
hold on; grid on;
plot(time,E,'-','LineWidth',1.6,'Color',[0 0.7 0]);
xlabel('\bfTime [sec]');
ylabel('\bfEnergy [J]');

linkaxes(ah,'x');
drawnow;

%% Animate result in a new plot
fh_anim = 2;
figure(fh_anim);
hold on; grid on;

% parameters of the pulley and cable wrapped around it (wrap length not
% accounted for in simulation!)
r_pulley = 0.1;
theta = 0:0.01:pi;
x_circ = r_pulley*cos(theta);
y_circ = r_pulley*sin(theta);

% plot each element, and get handles so positions can be adjusted in
% animation
plot(0,0,'.','MarkerSize',50,'Color',[0 0 0]);
ph_left = plot(-r_pulley*ones(1,2), nan(1,2),'-','LineWidth',3,'Color',[0.8 0 0]);
ph_arc = plot(nan,nan,'-','LineWidth',3,'Color',[0.8 0 0]);
ph_right = plot(r_pulley*ones(1,2), nan(1,2),'-','LineWidth',3,'Color',[0.8 0 0]);

% finish formatting axes
axis equal;
xlabel('\bfX');
ylabel('\bfY');
xlim([-0.4 0.4]);
ylim([-2.5 2*r_pulley]);

% animate each frame of results
saveFrameIdx = 0;
for tIdx = 1:anim_step:size(data,2)
    
    % extract state at current timestep
    s = data(1,tIdx);
    s_dot = data(2,tIdx);
   
    % update elements in current plot
    if(s <= sysParams.L/2)
        ph_left.YData = [s-sysParams.L/2 0];
        ph_arc.XData = x_circ;
        ph_arc.YData = y_circ;
        ph_right.YData = [0 -(sysParams.L/2+s)];
    else
        ph_left.YData = nan(1,2);
        ph_arc.XData = nan;
        ph_arc.YData = nan;
        ph_right.YData = (-s+(sysParams.L/2))+[0 -sysParams.L];
    end
    
    % update plot
    figure(fh_anim); % this is silly, needed to work around an issue with plot losing focus in MATLAB Online
    title(sprintf('Time: %6.3fs',time(tIdx)));
 	drawnow;
    pause(0.75*dt*anim_step);
    
end



% propagate state
function Xdot = stateProp(t,X,sysParams)

% recover parameters
L = sysParams.L;
g = sysParams.g;

% deconstruct state vector
s     = X(1);
s_dot = X(2);

% construct Xdot from differential equation
% note:     X    = [y y_dot]
% therefore Xdot = [y_dot y_ddot]
Xdot = zeros(2,1);
Xdot(1,:) = s_dot;

% equation of motion changes based on whether
% cable is unwinding around pulley or falling freely
if(s <= ???)            % TODO: REPLACE ??? WITH APPROPRIATE CONDITION
    Xdot(2,:) = ???;    % TODO: REPLACE ??? WITH EXPRESSION FOR s_ddot
else
    Xdot(2,:) =  ???;   % TODO: REPLACE ??? WITH EXPRESSION FOR s_ddot
end
end



