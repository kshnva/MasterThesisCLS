function S=ElectricEM(S)
% calculates the electrical source potential in the source elements of the
% electric network as a function of the shear stress in the corresponding
% hemodynamic network
% The source potential will be labeled sourceP in order to allow using 
% solvehemodyn to also solve this linear network


% a very rough first model
%for j=1:S.nESE
for j=1:S.nESEmiddle % limit this here to the mid-segment electrical sources
	WSS=S.IE (S.ESE(j).Hie).WSS;
    if S.PotassiumModel
        gK=S.ksens*abs(WSS)+S.gK0;
        gM=S.gNa+S.gLeak+gK; % total conductivity
        r= S.IE(S.ESE(j).Hie).r; % radius of the haemodynamic segment
        L=S.IE(S.ESE(j).Hie).l; % length of the hemodynamic segment
        S.ESE(j).Gs=4*pi*r*L*gM; % in paper this is Gm
        S.ESE(j).Ps=(S.gNa*S.ENa+S.gLeak*S.ELeak+gK*S.EK)./gM; % in paper this is Em
    else
    	S.ESE(j).Ps=S.E0-S.kMembraneShearSensitivity*(abs(WSS)-1); % hyperpolarization if WSS>WSSref
    end
end

% Flow-regulated capillary potentials (FullCouplingFlowRegulated model).
% Loops over the capillary sources only (indices nESEmiddle+1 : nESE).
if S.DoCapillaryCurrent && S.FlowRegulatedCapillaries
    for j = S.nESEmiddle+1 : S.nESE
        se = S.ESE(j).Hse;
        Q  = abs(S.SE(se).Qs);
        flowDeficit = max((S.QrefCapillary - Q) / S.QrefCapillary, 0);
        gKQ = S.kSensQ * flowDeficit;                              % [S/m2] flow-dependent K+ conductivity
        gsm = S.gCapBase + gKQ;                                    % [S/m2] total capillary membrane conductivity
        Esm = (S.gCapBase * S.ECapBase + gKQ * S.EK) / gsm;       % [V]    capillary equilibrium potential
        Gsm = S.kSmCapillary * gsm;                                % [S]    capillary membrane conductance
        % Thevenin equivalent at node: Gsc and Gsm in series
        S.ESE(j).Ps = Esm;
        S.ESE(j).Gs = S.GscCapillary * Gsm / (S.GscCapillary + Gsm);
    end
end