function S=defaultPars(topology,Sin)
% usage:  S=defaultPars(topology,[Sin]) 

% This defines the default model. The invar 'topology' indicates the
% topology to be set up

% Sin contains fields that override the default conditions, used for
% parameter variations 
% 

% All default parameters, modeling choices, specific values should go into
% this function. Changing the model from 
% default is only done in the GoProject main script. 

if nargin<2, Sin=struct; end

S.timestamp=datetime("now"); % keep for now, makes more sense to have this in the actual simulations, 
%is used in a rmfield statement when testing if the simulation already has been done

%% ===== 1) topology-independent parameters
S.model='24 October 2025'; % change this if the equations change!
S.annot='' ; % annotation of simulation

% === basic model parameters 
S.viscosity=4e-3;				% [Ns/m2]

% the model ODE expresses rdot as a function of the deviation of membrane
% potential from a reference value, originally rrdot=rr.*S.K3.*(S.Vref-Vm);
% now replace by rdot=r.*S.kReg.*(S.Vref-Vm)

S.tend=1e5;					% [s] 
%S.rstart=30e-6;					% [m] initial radii for searching solutions WSS==WSSref
S.r0=100e-6;% 30e-6						% [m] initial radii for starting ODE simulation

% === alternative models
S.AdaptationModel='WSSDrive'; % classes of models that use WSS as the drive for adaptation
%S.AdaptationModel='GeneralDrive'; % classes of models that use r^alpha * (deltaP/L)^beta as the drive for adaptation
S.AdaptAlpha=1; % alpha in adaptation rule, r^alpha in stimulus
S.AdaptBeta=1; % beta in adaptation rule (deltaP/L)^beta in stimulus
S.AdaptRateUp=1; % ratio of adaptation rate outward versus inward. Inward is kept constant 
S.scale_WSSref=1;				%	for heterogenous WSSref, use row vectors if not default

% === dynamics
%S.dynpar=""; % name of dynamic parameter
%S.freq=0; % frequency of possible oscillations of a parameter
%S.ampl=0; % amplitude
%S.phase=0; % phase in radials

% === obsolete, was used to find (unstable) initial equilibria
S.rSimplex=NaN;					%	radii resulting from fminsearch on WSS=WSSref, will be filled in during execution, if needed  [m]
S.InitDeviation=1.01;			%	relative radius at start of simulation [-]


%% adaptation model
S.kReg=0.01;					% [1/sV]
S.Vref=-0.05;					% [V], e.g. -50 mV as membrane potential that does not cause a change of radius. 
S.gainRadius2Vref=0; %200;			% [V/m] renders Vref dependent on radius, notes Nov 12 2025
S.offsetRadius2Vref= 50e-6;		% [m]


%% electrical communication parameters

% do either membrane source potential of potassium channel model
 S.PotassiumModel=1;
 S.DoElectricCoupling=1; 

% % % % == parameters associated with the source potential model - we are not
% % % % using this 
% % % % this is the model until 10/2025
% % % % Rc will be the coupling resistance, equaling rho*L/(2pi r * hc)= kElectricCoupling*L/r 
% % % % Gc would then be the inverse
% % % S.kElectricCoupling=1e8; %1			% [ohm]
% % % 
% % % % Rs is the membrane resistance (electrical source resistance) 
% % % % Rs=Rmembrane/(L*2pi r*2)= kElectricMembrane/ (L*r)
% % % S.kElectricMembrane=1;			% [ohm]
% % % 
% % % S.E0=0;							% [V], the membrane potential without shear stress and without currents to/from other segments
% % % S.kMembraneShearSensitivity=0.1;% [V/(N/m2)]


% == parameters associated with the potassium model and other changes in
% 11/2025
S.gNa=0.2; %3e-3;          % [S/m2]
S.ENa=60e-3;        % [Volt]
S.gLeak=0.2; %3e-3;       % [S/m2]
S.ELeak=0;          % [Volt] 
S.EK=-90e-3;        % [Volt]

S.ksens=0.06667; %1e-3;          % [S/m2 / N/m2] shear stress sensitivity of potassium channel opening
S.gK0=0.5; %0.75e-2;            % [S/m2] potassium channel conductance in absence of shear stress

S.ECheight=1e-6;    % [m] endothelial cell height
S.ECresistivity=0.1; %1e5;% [m/S] 

% == allow capillary currents
% this is the way to induce currents from the capillaries; for now the
% capillary conductance and potential are constant, this will need
% extension to a physiological model

S.DoCapillaryCurrent=0;			% boolean, include current generated in the capillaries or not
S.CapillaryPotential=-80e-3;	% [V]
S.CapillaryConductance=1;		% [S]

% == four coupling models (thesis): NoCoupling | ArtCoupling | FullCoupling | FullCouplingFlowRegulated
% DoElectricCoupling and DoCapillaryCurrent are set per-model in RunCouplingModels.
S.FlowRegulatedCapillaries=0;	% boolean: flow-dependent capillary K+ (FullCouplingFlowRegulated only)
% Capillary biophysical parameters (equations 9-12 in thesis)
S.gCapBase=0.2;					% [S/m²]    g_sb  baseline capillary membrane conductivity
S.ECapBase=-0.04;				% [V]       E_sb  baseline capillary equilibrium potential
S.kSensQ=0.2;					% [S/m²]    k_sensQ  capillary K+ sensitivity to normalised flow deficit
S.QrefCapillary=8e-14;			% [m³/s]    Q_ref reference capillary inflow (tune to network)
S.kSmCapillary=1e-6;			% [m²]      k_sm  area scaling: G_sm = kSmCapillary * g_sm
S.GscCapillary=2e-7;			% [S]       G_sc  fixed capillary coupling conductance

% === simulation parameters (as of 2024 02 09)
S.opts=[]; % ODE45 parameters, predefine the field to keep S fields consistent

% === parameters setting execution
S.KeepDynamics =1; % true if we want to keep the dynamics in state vars and derived vars, otherwise only the last timepoint is stored

% === parameters for making the videos of adaptation
S.videoFrameRate = 10; % Set frame rate
S.videoMaxLineThickness=5; % maximal line thickness in movie
S.videoMaxNFrames=inf; % maximum number of frames
S.videoQuality=100; % maximum quality
S.videoUsePNG=0; % set to one for decent or good quality, saving and reloading images, takes much time (1.5 s/frame versus 0.7 s/frame)






	
%% ===== 2) Definition of the topology
% 10-2025 several topologies have been removed for clarity, see older versions

% THE LEPTOMENINGEAL GEOMETRY
if strcmp(topology,'leptomeningeal_geometry')

S.ElementLength =2e-3; %1e-3;			% [m] 
S.devXY=0;						% [m] for honeycomb absolute deviation in position is devXY*(rand -0.5)  
S.ncombsx  =1;	%5 ;			% number of rows of hexagons, in x direction

% branchlevel is number of generations in the arterial trees perfusing the
% network; the terminal elements are connected to the leftmost and
% rightmost nodes respectively. 
S.branchlevel=4; % == 2025 

% S.ncomby is calculated from branchlevel below, because branchlevel may be
% overridden by non-default input parameters

% S.ncombsy  = 2^(S.branchlevel-2); % number of columns of hexagons, y direction
% branchlevel 2=> ncombsy=1; 3=>2, 4=>4 
% branchlevel 4=> ncombsy=4, 8 connections

% add randomness to the node positions
S.RandPosAmplitude=0;			% [m] 

% Number of penetrators per pial segment
S.npenetrator=1; %1; 

% depth of penetrators; they get a node at distal end so they are not SE
S.PenetratorDepth=5e-3;%3e-3;			% [m]

% define source and sink pressures
S.sourceP=80*133;				% [N/m2]
S.sinkP=15*133;					% [N/m2]

% 14-11-2025: define individual sink conductances,
S.GsSink=1.6e-14; % [m3/s] / [N/m2] conductance of a single sink

% Define the balance of total source and total sink conductance. Typically,
% this balance is in the order of 10. 10 means that the total source
% conductance is 10x higher than the total sink conductance.
S.balance_Gsin_Gsout=10;	% [-] balance of total input and output conductance

% all the individual source and sink conductances are derived from these
% two parameters and the number of sources and sinks. 

elseif strcmp(topology,'ThreeStar_geometry')
    S.sourceP=133*[100 70 20];
    S.sourceG=1e-12*[5 2 1]; 
else
	disp('unknown topology')
end


%% parameters that come with model variation
% these need to be defined here in order to keep S of each simulation the
% same structure

% dynamic parameters
S.dynpar=[];
S.ampl=[];
S.phase=[]; 
S.freq=[]; 


%% ===== handle non-default parameters for geometry, overriding the above defaults
% Loop through each field in Sin
fieldsSin = fieldnames(Sin);  % Get field names of Sin

for i = 1:numel(fieldsSin)
    field = fieldsSin{i};     % Get the field name
    
    % Replace field in S with field from Sin
    S.(field) = Sin.(field);
end

%% ==== any parameter that is calculated from other parameters needs to be defined here, since the base parameters may be non-default
if strcmp(topology,'leptomeningeal_geometry')
    S.ncombsy  = 2^(S.branchlevel-2); % number of columns of hexagons, y direction
    % branchlevel 2=> ncombsy=1; 3=>2, 4=>4 
    % branchlevel 4=> ncombsy=4, 8 connections
end

%% ===== 3a) Generation of the hemodynamic topology
S.topology=topology;
% fill the IE SE and IN fields as defined in the model
S=eval([S.topology '(S)']); % user defined connectivity model, can be replaced by other models

% determine the sizes for later quick reference
S.nIE=length(S.IE); % number of internal elements (segments), each connecting two nodes
S.nSE=length(S.SE); % number of connections to pressure sources/sinks from single node
S.nIN=length(S.IN); % number of nodes in the model, each node connects 2 or more internal or source elements

%% ===== 3b) Generation of the electrical topology
S=MakeElectricCircuit(S,0); % second parameter makes quick draw for test phase

%% ===== 4) Extension of initial values over the topology 

S.r0=S.r0*ones(1,S.nIE); % this is start set of radii

% obsolete, parameters for initial search for unstable equilibrium
%[S.IE.rstart]=deal(S.rstart); %*ones(1,S.nIE); % [-] initial relative radii for searching solutions WSS==WSSref
%S.rSimplex=S.rSimplex*ones(1,S.nIE); % radii resulting from fminsearch on WSS=WSSref
% keep this scalar, same for all elements
% S.InitDeviation = S.InitDeviation*ones(1,S.nIE); % initial deviation from
% unstable equilibria


