% restart
close all; clear; clc;

% parameters
sysParams.N_masses = 40;
sysParams.N_holdback = 20;   % # masses to hold at the origin until halfway through simulation
sysParams.m = 0.0005;       % [kg]
sysParams.k = 4;            % [N/m]
sysParams.c = 3.0;          % [Ns/m]
sysParams.g = 9.81;         % [m/s^2]

% simulation time parameters
t0 = 0;         % [s] simulation start time
tf = 8.0;       % [s] simulation end time
dt = 0.1;     % [s] timestep size

% initial conditions
x0 = zeros(2*sysParams.N_masses,1);
X = x0;

% data storage
time = t0;
data = x0;

% run simulation
for t = t0:dt:(tf-dt)
    
    % calculate timestep for ODE solving
    odeTime = [t t+dt];
    
    if(t > (tf/2))
        sysParams.N_holdback = 0;
    end

    % propagate state
    [T,X] = ode45(@(t,X) stateProp(t,X,sysParams),odeTime,X);
    X = X(end, :)';  % note: this step is necessary to keep state vector dimensions correct for next call to ode45()
    
    % store results from this timestep
    time(end+1)   = T(end);
    data(:,end+1) = X; % note: discarding state values at intermediate timesteps calculated by ode45()
end

%% animate results
figure;
hold on;
grid on; axis equal;
ph = plot(zeros(sysParams.N_masses+1,1),[0; nan(sysParams.N_masses,1)],'.','MarkerSize',20,'Color',[0.8 0 0]);
ylim([-1*max(data(sysParams.N_masses,:)),0]);
th = title(sprintf('%8.3fs',0));

for t_idx = 1:1:length(time)
    ph.YData = [0; -1*data(1:sysParams.N_masses,t_idx)];
    th.String = sprintf('%8.3fs',time(t_idx));
    drawnow;
    pause(0.1);
end

% propagate state
function Xdot = stateProp(t,X,sysParams)

N_masses = sysParams.N_masses;

% deconstruct state vector
% generate system matrix
T = toeplitz([-2 1 zeros(1,N_masses-2)]);
T(end,end) = -1;
A = [zeros(N_masses), eye(N_masses); (sysParams.k/sysParams.m)*T, (sysParams.c/sysParams.m)*T];
b = [zeros(N_masses,1); sysParams.g*ones(N_masses,1)];

A(1:sysParams.N_holdback,:) = 0;
A(N_masses + (1:sysParams.N_holdback),:) = 0;
b(N_masses + (1:sysParams.N_holdback),:) = 0;

Xdot = A*X+b;

end

