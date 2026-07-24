function [EIE,ESE]=ElectricConductance(S)
% calculates the electrical conductance of internal elements and source
% elements based on the radius and length, where radius is taken from the
% hemodynamic network

% okt 2025: replaced relative radii by absolute

% nov 2025: added the potassium model, where electrical conductance of the
% source already is calculated based on the potassium conductivity
EIE=S.EIE; ESE=S.ESE; % keep existing fields

if S.PotassiumModel
    % the source conductances based on the ionic conductances are
	% calculated in the ElectricEM function for the potassium model

	% calculate the coupling conductances
	if S.DoElectricCoupling
	   	for j=1:S.nEIE  
	    	r=S.IE (S.EIE(j).Hie).r; % get radius from the hemod network
        	l=S.EIE(j).l; % the electric model length, i.e. half the hemodynamic model length
        	EIE(j).G=2*pi*r*S.ECheight/(l*S.ECresistivity); 
    	end
	else
		[EIE.G]=vout(zeros(1,S.nEIE));
	end

else
    % below the old code





% copied from testelCoupling
% % % rl=S.rhoL/(2*pi*r*S.hEC); % long resistance per unit vessel length [ohm/m]
% % % RcHem=rl*L; % total long elect resistance of the hemodyn segment, [ohm]
% % % Rc=RcHem/2; % longitudinal resistance of the electrical segment [ohm], there are two such segments
% % % % we need G rather than 
% % % Gc=1/Rc; % G of an electrical internal element, there are two such elements that split this up. 
% % % 

% here the length is that of the electrical segment, so factor 2 already
% covered

% Calculation of coupling conductance, i.e. conductance in EIE, Rc in notes
for j=1:S.nEIE  
	r=S.IE (S.EIE(j).Hie).r; % get radius from the hemod network
    l=S.EIE(j).l;
% % % 	rl=S.rhoL/(2*pi*r*S.hEC);
% % % 	EIE(j).G=1/(rl*S.EIE(j).l); % electrical conductance
    EIE(j).G=r/(l*S.kElectricCoupling); 

end
% % % % source resistance, membrane resistance of the EC
% % % rm=0.5*S.Rm/(2*pi*r); % [ohm m], resistance for unit vessel length
% % % Rs=rm/L; % [ohm] total membrane resistance over the hemodyn vessel length, lower if longer vessel
% % % Gs=1/Rs; % there is only a single source element representing the full hemodyn vessel length
% % % 

% calculation of the membrane conductance, i.e. conductance of the ESE, Rs
% in notes
% the electrical source element has no length yet, use the length of the
% corresponding hemodyn element
%for j=1:S.nESE
for j=1:S.nESEmiddle % febr 2025 address only the hemodynamic sources in middle of segment, there are electric sources from the capillaire
	r=S.IE (S.ESE(j).Hie).r;
% % % 	rm=0.5*S.Rm/(2*pi*r);
% % % 	% we need to scale with the length of the original hemodyn segment!
	lhemod=S.IE(S.ESE(j).Hie).l;
	ESE(j).Gs=lhemod*r/S.kElectricMembrane; % use rGs for now, otherwise not recognized in solvehemodyn2024
end

end 
