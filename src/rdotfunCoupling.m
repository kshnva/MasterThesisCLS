function rdot=rdotfunCoupling (t,r,S)
% calculates the time derivative of the state vars (the vector of relative
% radii).
% t: time (s), used for the dynamic load model
% r: vector of radii (m)
% par: structure with the model parameters
% rdot: vector of time derivatives (m/s)
% this is the coupled model, where WSS drives membrane potential
% and membrane potential drives rrdot

%% 2024 S.IE structure: deal r over S.IE.r
[S.IE.r]=vout(r);


%% DYNAMIC LOAD: modify a parameter if one is chosen to be dynamic
% dynpar can be combined with any model, it sets the parameters before
% applying the adaptation model
%if isfield(par,'dynpar')
if ~isempty(S.dynpar)
   	instruction=['par.' S.dynpar '=par.' S.dynpar '.*(1+par.ampl.*sin(par.phase+2*pi*par.freq*t));'];
	eval(instruction);
    %this won't work because pars may be vectors, rPs(2) is not a field
	%name:
	%par.(par.dynpar)=par.(par.dynpar).*(1+par.ampl.*sin(par.phase+2*pi*par.freq*t));
end

%% WSS drives membrane potential drives adaptation 2025 absolute values
% step 1: calculate the wall shear stresses

[S]=calcConductance2025(S);
[S]=solvehemodyn2025(S);
[S]=calcWSS2025(S);

% step 2: calculate the equilibrium potentials based on the shear
% stresses

% 11-2025 this includes calculating the membrane conductance, if we use the
% potassium conductance model
[S]= ElectricEM(S); 

% step 3: solve the electrical network, deriving V in the middle of the
% segments, which is the drive for adaptation. We use solvehemodynamics
% for this, but this is actually applied to the electrical system

[EIE,ESE]=ElectricConductance(S); 
S.EIE=EIE; S.ESE=ESE; 
[S]=solvehemodyn2025(S,'electrical');

%env=[S.ESE.node]; % these are the connected electrical nodes, halfway the hemodynamic elements
% febr 2025: only do this on electrical sources halfway the segments, not
% the capillary sources
env=[S.ESE(1:S.nESEmiddle).node]; % these are the connected electrical nodes, halfway the hemodynamic elements
Vm=[S.EIN(env).P]'; % these are the actual membrane potentials 

% step 4: calculate the temporal derivatives of the radii
%rdot=r.*S.kReg.*(S.Vref-Vm);
rdot=r.*S.kReg.*(S.Vref-S.gainRadius2Vref*(r-S.offsetRadius2Vref)-Vm);

%% DIRECTION SENSITIVITY
% direction sensitivity can be combined with any model, it scales rdot
% after applying the adaptation model
if S.AdaptRateUp~=1
    rdot=((S.AdaptRateUp-1).*(rdot>0)+1).*rdot;
end	

%% HYSTERESIS obsolete
% % % % hysteresis can be combined with all WSS models
% % % if isfield(par,'HystSpan')
% % % 		DoChange=(rAdaptationDrive<1-par.HystSpan) | (rAdaptationDrive>1+par.HystSpan); % only have a change in radius if rWSS is far enough from unity
% % % 		%rrdot=DoChange.*rr.*(rWSS-1) * par.K3;
% % %         rrdot=DoChange.*rrdot;
% % % end 
end



