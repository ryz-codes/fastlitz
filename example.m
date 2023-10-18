%% Generates the simulation model via command-line 
% Command-line is more powerful, but the GUI is easier to use.
% Everything in this section can be done in the GUI, by clicking
% "Initialize Sim. Model".

% Definitions
% Syntax for the litz.wire command:
% > geom = litz.wire(len,pitch,ns,srad,insu,drad,odd)
len = 0.01; %[m] length of wire
pitch = [1e-2, 2e-2]; %[m]
ns = [3, 3]; % number of strands at each level
srad = 1e-4/3; %[m] strand radius
insu = 0.1; % (diam + insulation) / diam - 1
drad = 3; % discretization param
odd = true; % filament at center?

% Generate model
wire_model = litz.wire(len,pitch,ns,srad,insu,drad,odd);

% Render the model
litz.showGeom(wire_model);

% Set up the simulation model. This takes a few minutes.
model = pmag(wire_model);

%% Run a variety of simulations using the generated model
% Run sweeps. We can run as many sweeps as we like while avoiding repeating
% the setup.
f = logspace(3,7,20); % [Hz]

% To extract terminal impedance, use
%   [Z] = model.impedance(f)
Z = model.impedance(f);

% To compute power loss due to an incident magnetic field, use
%   [P] = model.power(f,[Bx,By,Bz]) % incident B field
%   [P] = model.power(f,[Hx,Hy,Hz],'A/m') % incident H field
P = model.power(f,[0,0,1],'A/m');

% Plot
figure;
loglog(f,real(Z));
figure;
loglog(f,P);