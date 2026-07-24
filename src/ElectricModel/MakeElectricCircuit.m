function S=MakeElectricCircuit(S,MakeQuickDraw)
% usage S=MakeElectricCircuit(S,[MakeQuickDraw])
% this functions adds an electric circuit based on the existing structure
% of internal nodes and internal elements.
% each IE is represented by two Electrical Internal Elements and an
% Electrical Internal Node
% Each EIN connects the two EIE to an ESE (electrical source element), of
% which the driving voltage and source conductance can be based on the
% shear stress.
% This function defines the connectivity and length, as well as the
% references between both trees. E and Gs are set at initialization and
% during the adaptation process
if nargin==1, MakeQuickDraw=0; end

ieie=0; % counter for electrical internal elements
iese=0; % counter for electrical source elements
iein=S.nIN; % counter for electrical internal nodes 
for ie=1:S.nIE
	iein=iein+1; % new electrical internal node index
	ieie=ieie+1; % new electrical internal element, representing first half of hemodynamic element
	S.EIE(ieie).nodes(1)=S.IE(ie).nodes(1); % start from same node number as hemodynamic model
	S.EIE(ieie).nodes(2)=iein; % reference to new node
	ieie=ieie+1; % new electrical internal element, second half
	S.EIE(ieie).nodes(1)=iein; % reference to new node
	S.EIE(ieie).nodes(2)=S.IE(ie).nodes(2); % connect end of this element to end of hydrodynamic model
	[S.EIE((ieie-1):ieie).l]=deal(S.IE(ie).l/2); % set the lengths
	
	iese=iese+1; % new source element
	S.ESE(iese).node=iein;
	
	% refer to this electrical elements and nodes in the hemodyn network
	S.IE(ie).Eie=[ieie-1, ieie]; 
	S.IE(ie).Ein=iein; % the middle internal nodes, the other two are in nodes anyway
	S.IE(ie).Ese=iese; % the electrical source element, needed because its E will depend on shear stress in IE
	
	% and in the electrical source and internal elements refer back to the
	% hemodynamic one
	[S.EIE((ieie-1):ieie).Hie]=deal(ie);
	S.ESE(iese).Hie=ie;
end

% 28 Febr 2025: on top of the above 'middle of vessel' shear sensitive
% electrical source elements, we want to have external potentials
% representing the capillaries (the sinks). Continue building ESE's but
% remember which of these are 'middle of vessel' 
S.nESEmiddle=iese;

if S.DoCapillaryCurrent
for se=1:S.nSE
	if ~S.sources(se)% select the sinks and not sources
		iese=iese+1; % new electrical source element at this hemodynamic source element. 
		node=S.SE(se).node; % the hydraulic node
		S.ESE(iese).node=node; % new electrical source element connected to this electrical node
		S.ESE(iese).Ps=S.CapillaryPotential;
		S.ESE(iese).Gs=S.CapillaryConductance; % for now scalar constants
		S.ESE(iese).Hse=se; % back-ref: hemodynamic sink SE index (used by ElectricEM for flow regulation)
		% % % S.EIN(node).nsources=S.EIN(node).nsources+1; % dit zit in
		% MakeNodeTable
		% % % S.EIN(node).se=iese; 
	end
end
end

% have the counters for future reference
S.nEIE=ieie;S.nESE=iese;

% make the electrical node table
[S.EIN,S.nEIN]=MakeNodeTable(S.EIE,S.ESE); %% maybe relable the ie to eie now


% define the positions of the electrical nodes

for ein=1:S.nIN
	S.EIN(ein).pos=S.IN(ein).pos; % these were the identical-numbered nodes at the ends
end

for ein=S.nIN+1:S.nEIN
	% this node should be halfway the connected nodes, which are both
	% identical to the hemodynamic nodes so we know their position
	%ie=S.EIN(ein).ie; % the connected internal elements 
	S.EIN(ein).pos=0.5*S.IN(S.EIN(ein).cn(1)).pos+0.5*S.IN(S.EIN(ein).cn(2)).pos;
end

%somehow we use pos and (x,y) combined, needs to be sorted out, we need x y for plotting now
for ein=1:S.nEIN
	S.EIN(ein).x=S.EIN(ein).pos(1);S.EIN(ein).y=S.EIN(ein).pos(2);
end


%% quick draw to see what we get
if MakeQuickDraw
    hf=figure(555); clf;hold on; title('quick draw ELECTRIC configuration in honeycomb geometry definition')
    for i=1:S.nEIN
	    plot(S.EIN(i).x,S.EIN(i).y,'og')
    end
    for j=1:S.nEIE
	     xx=[S.EIN(S.EIE(j).nodes).x];yy=[S.EIN(S.EIE(j).nodes).y]; %% xx are x values of left and right node
	     plot(xx,yy,'-b')
    end
    v=find([S.EIN.nsources]); % these indices in IN are connected to SE
    for k=1:length(v)
	     plot(S.EIN(v(k)).x,S.EIN(v(k)).y,'or')
    end
    hold off
end

