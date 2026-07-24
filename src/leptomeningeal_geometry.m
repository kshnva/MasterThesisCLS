function S = leptomeningeal_geometry(Sin)
% usage: S = leptomeningeal_geometry(Sin)
% Generates a leptomeningeal circulation, with two entrances, a pial
% collateral network and penetrators

% keep the current fields in the model
S=Sin;

% we will loose segments, having their conductivity towards zero and
% generating matrix warning. Disable these warnings but remember we did so
% mmm... I still get the warnings
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
S.SingularMatrixWarning='off';

%% Define the version of this geometry
% note that this file contains choices and parameters, if these change and if these changes do not appear
% in S, a new version number is needed.
% new version number is need.
S.geom.version='leptomeningeal version October 23 2025';
%% === THE CASE ==============
% leptomeningeal geometry, entrance on left could be from MCA and on right from ACA
% branching trees from both sides that are connected to a honeycomb
% geometry. This whole thing is the pial circulation. 
% Each pial vessel then has one or more penetrators, each consisting of an
% IE in series with a sink element (SE). 

% % % % if we haven't defined this yet, define the size of the network
% % % if ~isfield(S,'ncombsx')
% % % 	S.ncombsx  = 1; % number of collateral hexagonal loops in x direction 
% % % 	S.branchlevel=2; % number of generations MCA and ACA tree, 2=> 3 segments, one mother and two daughters
% % % 	S.ncombsy  = 2^(S.branchlevel-2); % same in y direction
% % % 	% branchlevel 2=> ncombsy=1; 3=>2, 4=>4 
% % % end

%% build a left and right tree with sources
% the inputs are maintained, sinks are later cut off when fused to the
% collateral plexus
S2=S; 
S2=tree_geometry(S,1); % second argument: use this geometry as subgeometry 
S3=S;
S3=tree_geometry(S,1);
% % % % exclude the geometry of the source elements, otherwise the solver tries
% % % % to calculate source element WSR
% % % S2.SE=rmfield(S2.SE,{'r','l'}); S3.SE=rmfield(S3.SE,{'r','l'});
for in=1:S2.nIN, S2.IN(in).pos=[0; 0]; S2.IN(in).side='L';end
for in=1:S3.nIN, S3.IN(in).pos=[0; 0]; S3.IN(in).side='R';end

%% build the honeycomb
S=honeycomb_geometry(S,1); % second argument: use this geometry as subgeometry 
% we won't be needing the source and sink connections
S.SE=S.SE([]);
S.nSE=0;

% % % try % SE had these fields in 2022 but not 2025
% % %     % also remove the l and r fields, they force calculation of SE wall shear
% % %     % stress in the solver. 
% % %     S.SE=rmfield(S.SE,{'r','l'}); 
% % % end


% the node table is reconstructed later, no need to adapt it

%% connect the trees to the honeycomb
% left (note that nodes are renumbered after first tree is connected!)

% figure out how many nodes need to be connected
NFirstSecondColumn=[4 5 8 9 12 13 16 17 20 21]; % these number of nodes are in first and second column for ncombsy = 1,2,3 etc *** make this into rule
NFSC=NFirstSecondColumn(S.ncombsy); % for this specific network size

% figure out which nodes on the left side should be connected
nodes_position=[S.IN.pos];
nodes_x=nodes_position(1,:);
[~,I]=sort(nodes_x); 
Connect_Left=I(1:NFSC);  % possible connection points: the two leftmost columns in the honeycomb

% exclude the extreme y position(s), ncombsy odd: exclude at both ends,
% ncombsy even: exclude lowest y position
nodes_y=nodes_position(2,Connect_Left);%[S.IN(Connect_Left).y];
[~,J]=sort(nodes_y);
Connect_Left=Connect_Left(J); % now sorted on Y position
if mod(S.ncombsy,2) % uneven
	Connect_Left=Connect_Left(2:NFSC-1);
else
	Connect_Left=Connect_Left(2:NFSC);
end

% join the left tree and the honeycomb
S=jointrees(S,Connect_Left,S2,[S2.SE(2:end).node]); % nodes are removed from the second network, the first network already has coordinates of the node

% figure out which nodes on the right side should be connected
% note that the nodes have been renumbered above
[~,I]=sort(nodes_x); 
Connect_Right=I(end-NFSC+1:end);  % possible connection points: two most right columns

% exclude the extreme y position(s), ncombsy odd: exclude at both ends,
% ncombsy even: exclude highest y position
nodes_y=nodes_position(2,Connect_Right);%[S.IN(Connect_Left).y];
[~,J]=sort(nodes_y);
Connect_Right=Connect_Right(J); % now sorted on Y position
if mod(S.ncombsy,2) % uneven
	Connect_Right=Connect_Right(2:NFSC-1);
else
	Connect_Right=Connect_Right(2:NFSC);
end

% join the right tree as well 
S=jointrees(S,Connect_Right,S3,[S3.SE(2:end).node]); 

%% match the coordinates of the tree nodes to those of the honeycomb 
% The y coordinate of each mother is the average of that of her daughters
% The x coordinate: alike, but shifted to left for the left tree, and to
% the right for the right tree. 
prox_nodes=[S.IE.nodes];prox_nodes=prox_nodes(1:2:end); % proximal nodes of each IE
dist_nodes=[S.IE.nodes];dist_nodes=dist_nodes(2:2:end); % distal nodes of each IE (i.. closer to the honeycomb)
% only ie belonging to the trees have branch level>0; work from distal to
% proximal, simultaneously for both trees in this loop
for gen=S.branchlevel:-1:1
	tree_prox_nodes=prox_nodes(find([S.IE.generation]==gen));
	tree_dist_nodes=dist_nodes(find([S.IE.generation]==gen));
	[S.IN(tree_prox_nodes).ndist_ie]=deal(0); 
	for j=1:length(tree_prox_nodes) % add all the positions of distal nodes to the proximal node
		if strcmp(S.IN(tree_prox_nodes(j)).side,'L'), k=-1; else k=1; end % horizontal shift of mother versus daughters, negative for left tree
		S.IN(tree_prox_nodes(j)).pos(1)=S.IN(tree_prox_nodes(j)).pos(1)+S.IN(tree_dist_nodes(j)).pos(1)+k*S.ElementLength*2^0.5; % this works towards average of both distal node x values, shifted by the width of a hexagonal
		S.IN(tree_prox_nodes(j)).pos(2)=S.IN(tree_prox_nodes(j)).pos(2)+S.IN(tree_dist_nodes(j)).pos(2); % this works towards average of both distal node y values
		S.IN(tree_prox_nodes(j)).ndist_ie=S.IN(tree_prox_nodes(j)).ndist_ie+1; % number of distal elements to this node
	end
	unique_tree_prox_nodes=unique(tree_prox_nodes);
	for j=1:length(unique_tree_prox_nodes) % the above loop summed the daughter positions, now divide by the number of daughters
		S.IN(unique_tree_prox_nodes(j)).pos=S.IN(unique_tree_prox_nodes(j)).pos/S.IN(unique_tree_prox_nodes(j)).ndist_ie; % normalize the proximal node positions
	end
end

%% add randomness to the x and y positions
xrand=(rand(1,S.nIN)-0.5)*S.RandPosAmplitude*S.ElementLength;
yrand=(rand(1,S.nIN)-0.5)*S.RandPosAmplitude*S.ElementLength;
for in=1:S.nIN, S.IN(in).pos(1)=S.IN(in).pos(1) +xrand(in); end
for in=1:S.nIN, S.IN(in).pos(2)=S.IN(in).pos(2) +yrand(in); end

%% define these IE as pial
[S.IE.pial]=deal(1); 

%% define the third dimension to allow for the penetrators towards the z direction
for in=1:S.nIN, S.IN(in).pos=[S.IN(in).pos; 0]; end

%% introduce the penetrators
% each of the pial vessels, including the left and right supplying trees,
% will get penetrators. S.npenetrator should be defined in default
% hf=QuickDraw(S); use this to see what we have before we start adding
% penetrators
hnIE=S.nIE; % current number of internal elements
for ie=1:hnIE
	nodes=S.IE(ie).nodes;
	dpos=(S.IN(nodes(2)).pos-S.IN(nodes(1)).pos)/(S.npenetrator+1); % shift in x y for penetrator relative to first node
	for pe=1:S.npenetrator
		%S=ConnectPenetrator(S,ie,dpos*pe);
		% splitting occurs from right to left
		S=ConnectPenetrator(S,ie,dpos*(S.npenetrator-pe+1));
	end
end

%% make the node table
% above the node numbers are defined in IE, now adapt the IN structure
% related to connectivity: nconnect nsources cn ie and se 
% cn vector of connected nodes, ie same size vector of ie numbers forming
% the connection, nsources and se the connected source elements
[S.IN,S.nIN]=MakeNodeTable(S.IE,S.SE,S.IN); 

%% Set the element lengths
% we moved and split the elements, now set their lengths such that the
% conductances are calculated correctly
S.IE=LengthFromPosition(S.IE,S.IN);

%% Define the source elements
% these elements are defined above but these need to get their pressures
% and conductances

% % % % the inlet to the two trees
% % % 
% % % Gin=S.Gstot*S.sourceG_sinkG_balance/2; % total source conductance equally distributed
% % % [S.SE(1:2).Gs]=deal(Gin); 
% % % 
% % % % the outlets distal of the penetrators
% % % 
% % % Gout=S.Gstot*(1-S.sourceG_sinkG_balance)/(S.nSE-2); % total sink conductance equally distributed
% % % [S.SE(3:end).Gs]=deal(Gout);



% 14-11-2025
n_in=2;
n_out=S.nSE-n_in;

[S.SE(1:n_in).Ps]=deal(S.sourceP); % S.SourceP is the high pressure, the others get low pressure
[S.SE(n_in+1:end).Ps]=deal(S.sinkP);

Gsinktot=S.GsSink*n_out; % total sink conductance
Gsourcetot=S.balance_Gsin_Gsout*Gsinktot; % total source conductance
[S.SE(1:n_in).Gs]=deal(Gsourcetot/n_in); 
[S.SE(n_in+1:end).Gs]=deal(S.GsSink);

S.sources=zeros(1,S.nSE); % these are used by calcXdot, also to allow steps in input pressure
S.sources(1:n_in)=1; 

%% define the default initial state of the model
% the definition and number of state vars in each IE depend on the model (e.g. 1-state, 5-state)
% so this part of the code is different for different models

% if strcmp(S.model(1:3),'One')
% 	% initial state definition for 1-state models: radius as defined in the
% 	% defaultOneStatePars model 
% 	S.X0=ones(length(IE),3)* S.r0;		% radii (m)
% 		
% % % 	still have to figure out whether or not have r here in network, or
% % have X0 distributed over network @@@@@@
% 
% elseif strcmp(S.model, 'FiveState version sept 16 2021')
% 	% initial state definition for 5-state models    
% %       S.X0 = [S.r0 * ones(length(IE),1);  S.tone0 * ones(length(IE),1);  S.span0 * ones(length(IE),1);...
% %       S.rmslack0 * ones(length(IE),1); S.wCSA0 * ones(length(IE),1)]'; 
% %      S.X0 = [S.strain0 * ones(length(IE),1);  S.tone0 * ones(length(IE),1);  S.span0 * ones(length(IE),1);...
% %      S.rmslack0 * ones(length(IE),1); S.wCSA0 * ones(length(IE),1)]'; 
% 
%     %S.X0 = [S.r0;  S.tone0;  S.span0; S.rmslack0; S.wCSA0];
% end
%% remove temporary fields from the nodes
S.IN=rmfield(S.IN,{'side','ndist_ie'});
%% quick draw to see what we get
% hf=QuickDraw(S,234,'quick draw of configuration in leptomen geometry definition');

%% SUB FUNCTIONS

function S=ConnectPenetrator(S,ie,dpos)
	% splits up a pial vessel in two pial vessels and connects a penetrator
	% and a sink element
	
	% split pial vessel
	hn=S.IE(ie).nodes(2);
	S.IE(ie).nodes=[S.IE(ie).nodes(1) S.nIN+1]; % new right node original ie
	S.IE(S.nIE+1)=S.IE(ie); % new ie, with same field values as ie for now
	S.IE(S.nIE+1).nodes=[S.nIN+1 hn ]; %new nodes new ie
	
	% define the new node position at origin of penetrator
	% note we split the segment coming from the first node over and over
	ln=S.IE(ie).nodes(1);
	rn=S.IE(S.nIE+1).nodes(2); 
	S.IN(S.nIN+1)=S.IN(S.nIN); % make a new node, we need its definition but all the connectivity comes in MakeNodeTable
	%S.IN(S.nIN+1).pos=0.5*(S.IN(ln).pos+S.IN(rn).pos);
	S.IN(S.nIN+1).pos=S.IN(ln).pos+dpos;
	%S.IN(S.nIN+1).pos=S.IN(rn).pos-dpos;	
	
	% new penetrator
	S.IE(S.nIE+2)=S.IE(ie); % new ie, with same field values as ie for now
	S.IE(S.nIE+2).nodes=[S.nIN+1 S.nIN+2]; % new distal node
	S.IE(S.nIE+2).pial=0; % meanINg it is a penetrator
	
	% set the properties of this penetrator
	% ==> to do, for now
		%S.IE(S.nIE+2).diam=30e-6; % to do, diameters are not defined at all yet **
	
	% define the new node position at distal end of penetrator
	S.IN(S.nIN+2)=S.IN(S.nIN); % make a new node
	S.IN(S.nIN+2).pos=S.IN(S.nIN+1).pos + [0; 0; -S.PenetratorDepth]; % same x-y and negative z, into tissue
	% new sink element
	S.SE(S.nSE+1).node=S.nIN+2; 
	% still have to do other SE fields
		
	% update number of nodes and elements
	S.nIN=S.nIN+2;
	S.nIE=S.nIE+2;
	S.nSE=S.nSE+1;
end




end
