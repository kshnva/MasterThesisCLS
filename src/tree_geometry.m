function [S] = tree_geometry(Sin,subgeometry)
% generates a tree geometry (for now simple symmetric geometry, same number
% of generations along all paths)
% Sin.branchlevel indicates number of generations,  1=> 1 segments, 2=> 3 segments, one mother and two daughters
% if subgeometry is defined and true, this function is used as part of a
% larger network generator. Several parts are then not executed
% see usage in leptomeningel_geometry

if nargin==1
	subgeometry=0;
end

% keep the current fields in the model
S=Sin;

if ~ subgeometry
% we will loose segments, having their conductivity towards zero and
% generating matrix warning. Disable these warnings but remember we did so
% mmm... I still get the warnings
%warning('off','MATLAB:singularMatrix')
%warning('off','MATLAB:nearlySingularMatrix')
%S.SingularMatrixWarning='off';
end 
%% Define the version of this geometry
if ~ subgeometry
% note that this file contains choices and parameters, if these change and if these changes do not appear
% in S, a new version number is needed.
% new version number is need.
S.geom.version='tree version 2025 10 23';
end
%% === THE CASE ==============
% if we haven't defined this yet, define the branch level
% % % if ~isfield(S,'branchlevel')
% % % 	S.branchlevel=3; % 4 end segments
% % % end


%% define Internal Elements

% generate the first segment
S.IE=[]; % otherwise the previous fields remain
S.IN=[];
S.SE=[];

S.nIE=1;
S.nIN=2; 
S.nSE=1;

S.IE(1).generation=1;
S.IE(1).nodes=[1 2];
S.SE(1).node=1;

S=AddSubtrees(S); % add left and right subtrees

	function S=AddSubtrees(S)
		mie=S.nIE;
		if S.IE(mie).generation<S.branchlevel
			% add two daughters
			for j=1:2 % for two daughters
				S.nIE=S.nIE+1; % generate a new daughter
				S.IE(S.nIE).nodes(1)=S.IE(mie).nodes(2); % define proximal node of daughter
				S.nIN=S.nIN+1;
				S.IE(S.nIE).nodes(2)=S.nIN; % new distal node
				S.IE(S.nIE).generation=S.IE(mie).generation+1; % one generation higher than mother
				S=AddSubtrees(S); % recursive addition of subtree
			end
		else
			% add sink to distal node of non-mother 
			S.nSE=S.nSE+1;
			S.SE(S.nSE).node=S.nIN;
		end
	end


%% generate the node table / don't change this
% we now have a table of elements; we also need to have a table of internal nodes 
% this function makes a table of nodes, each node has information on the
% connected segments and will have pressure data
% never change this!
[S.IN,S.nIN]=MakeNodeTable(S.IE,S.SE);
%% define the position of the nodes
% from left to right, x defined by generation of mother, y by number of already assigned positions, get these number into the IN structure

genCount=zeros(1:S.branchlevel);
S.IN(1).gen=0; % proximal of first IE, connects to SE
S.IN(1).genCount=1;
for in=2:length(S.IN)
	gen=S.IE(S.IN(in).ie(1)).generation;
	S.IN(in).gen=gen;
	genCount(gen)=genCount(gen)+1;
	S.IN(in).genCount=genCount(gen);
end
v=1:S.branchlevel+1;
xv=v.^0.5;
yv=1:2^(S.branchlevel);

% keep for now but these values are redefined if this tree is used as part
% of a larger geometry

for in=1:length(S.IN)
	S.IN(in).pos(1,1)=xv(S.IN(in).gen+1); %  Oct 2025: now a column vector
	S.IN(in).pos(2,1)=yv(S.IN(in).genCount);
end

%% have the separate x and y => 2025: obsolete
% % % for in=1:length(S.IN)
% % % 	S.IN(in).x=S.IN(in).pos(1);
% % % 	S.IN(in).y=S.IN(in).pos(2);
% % % end


%% define the lengths and other properties of source-connecting elements
% we need to define length and radius and calculate conductance
% or alternatively just enter conductance
%[SE.l]=deal(ElementLength); % L was defined above for determining positions, we will use it for the source element lengths too 
%[SE.r]=deal(S.r0); % all source elements have the same radius as the starting radius of the IE

[S.SE.Gs]=deal(0); % As part of the leptom model, this is not used. For stand-alone to be defined ***

%% define the lengths of the internal elements / don't change this
% calculate the lengths of the segments from the positions of the nodes,
% assuming straight segments
% might be obvious in some cases like the wheatstone, but not in others

% this is redefined in leptom geometry?
% % % gen=[S.IE.generation];
% % % L=gen.^-0.7.*(1+0.3*(rand(1,S.nIE)-0.5));sumL=sum(L);
% % % [S.IE.l]=vout(L/sumL);

% % % S.IE=LengthFromPosition(S.IE,S.IN); 
% % % [S.IE.l]=vout([S.IE.l]/sum([S.IE.l]));
%% define the external pressures
% note that these values might be changing in time, if they are chosen as
% input variables. In that case it should not be necessary to initialize
% them here. 

% this is redefined in leptom geometry?
% % % % the SE connected to the input node gets the high pressure, the others low
% % % % presures
% % % [S.SE(1).Ps]=deal(S.sourceP);	% (N/m2)
% % % [S.SE(2:end).Ps]=deal(S.sinkP);		% (N/m2)

%% define the default initial state of the model
% still done?

%% assign to the output structure
% nIEt nodig, S.nIN en zo heb ik ook al
%% quick draw to see what we get
%hf=QuickDraw(S,233,'quick draw of configuration in tree geometry definition');

% % % hf=figure(234); clf;hold on; title('quick draw of configuration in tree geometry definition')
% % % for i=1:S.nIN
% % % 	plot(S.IN(i).x,S.IN(i).y,'og')
% % % end
% % % for j=1:S.nIE
% % % 	 xx=[S.IN(S.IE(j).nodes).x];yy=[S.IN(S.IE(j).nodes).y]; %% xx are x values of left and right node
% % % 	 plot(xx,yy,'-b')
% % % end
% % % v=find([S.IN.nsources]); % these indices in IN are connected to SE
% % % for k=1:length(v)
% % % 	 plot(S.IN(v(k)).x,S.IN(v(k)).y,'or')
% % % end
%%
end