function [S] = honeycomb_geometry(Sin,subgeometry)
% usage: S=honeycomb(S,[subgeometry])
% Honecomb generates a honeycomb-like geometry 
% Further explanation is in the wheatstone_geometry file
% if subgeometry is defined and true, this function is used as part of a
% larger network generator. Several parts are then not executed
% see usage in leptomeningeal_geometry
% okt 2025: removed x and y fields of nodes, keep position as 2D column
% vector

if nargin==1
	subgeometry=0;
end

% keep the current fields in the model
S=Sin;

if ~ subgeometry
% we will loose segments, having their conductivity towards zero and
% generating matrix warning. Disable these warnings but remember we did so
% mmm... I still get the warnings
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
S.SingularMatrixWarning='off';
end 
%% Define the version of this geometry
if ~ subgeometry
% note that this file contains choices and parameters, if these change and if these changes do not appear
% in S, a new version number is needed.
% new version number is need.
S.geom.version='honeycomb version 2025 10 25';
end
%% === THE CASE ==============
% Honeycomb structure consisting of hexagonals, length and width as defined
% in the parameters hexx and hexy

% if we haven't defined this yet, define the size and number of inflow and
% outflow vessels. (assumed this is always defined if this is subgeometry)
if ~isfield(S,'ncombsx')
	S.ncombsx  = 3; 
	S.ncombsy  = 1;
end

% define derived size indicators to get the configuration right
hexx=2+2*S.ncombsx;
hexy=S.ncombsy; 


%% define Internal Elements
for i = 3:1:hexx
    IE(i-2).nodes = [i-2, i-1];   % 1~(hexx-1)
    if hexy ==1  % if only one hexagon in y direction
        IE(i-3+hexx-1).nodes = [i-3+hexx, i-3+hexx+1];
        for k = 0:1:(hexx/2)-1
            IE(hexx*hexy + hexx -3 + k).nodes = [2*k+1, hexx+(2*k)]; % vertical element
        end
    else % more than one hexagonal
        for n = 1:1:hexy % y direction, parallel element
            IE(i-3+(n*(hexx-1))).nodes=[i-3+(n*(hexx)),i-3+n*(hexx)+1]; % nth (hexx-1)~(n+1)th (hexx-1)
            if n ~= hexy %for shift element
                IE((hexx-2)*(n+1)+n).nodes = [(hexx-1)*(n+1)+n-1, (hexx-1)*(n+1)+n];
            end
        end
        for j = 0:1:hexx/2-1 % vertical element
            for k = 1:1:hexy
                if hexy ==2
                    IE(j+1+ 1/2*hexx*(k-1) +(hexx -2)*hexy + hexx-1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2];
                else
                    if k ==1
                        IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2];
                    elseif k ==2
                        IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2+1];
                    elseif mod(k,2)==0 %even
                        if k == hexy
                            IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2];
                        else
                            IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2+1];
                        end
                    else %odd
                        IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+j*2, hexx*k+j*2];
                    end
                end
            end
        end
    end
end

%% define node connection table part one - the internal elements
% we normally would use the statement [IN,nin]=MakeNodeTable(IE,SE);
% but we now first need the nodes before we hook on the source elements
% so this is the first part of MakeNodeTable
nIE=length(IE);

% total number of internal nodes
nIN=max([IE.nodes]);  

%no connections at start
for i=1:nIN 
    IN(i).nconnect=0;
end

% fill the internal connections based on the definition of IE,  
for ie=1:nIE
    node(1)=IE(ie).nodes(1); 
    node(2)=IE(ie).nodes(2);
    for j=1:2
        k=node(j);
        IN(k).nconnect=IN(k).nconnect+1;
        nc=IN(k).nconnect;
        IN(k).cn(nc)=node(3-j); % 3-j means 1=>2 and 2=>1
        IN(k).ie(nc)=ie; % points back to element table
    end
end


%% define the position of the nodes
% in this model, the hexagons are regular, all IE have the same length,
% taken unity here and scaled at end

% ED March 18, 2022: well..., they are not regular. I never realized this, but these are hexagons with
% angles 135-90-135-135-90-135 rather than equal angles of 120 
% This cannot be solved by strethching in x and shrinking in y, because the
% lengths then become unequal. => keep it this way

for i=0:1:((hexx)/2)-1              %i=0,1,2,3
    IN(1+2*i).x=2*i;                  %1,3,5,7
    IN(1+2*i).y=1;                    %1,3,5,7
    for j=0:1:((hexx)/2)-2          %j=0,1,2
        IN(2+2*j).x=1+2*j;            %2,4,6
        IN(2+2*j).y=0;                %2,4,6
        for k=0:1:hexy-1            %k=0,1
            if k==0
                IN(hexx+2*i).x=2*i;                         %8,10,12,14
                IN(hexx+2*i).y=1+sqrt(2);
                IN(hexx+1+2*j).x=1+2*j;                     %9,11,13
                IN(hexx+1+2*j).y=2+sqrt(2);
            elseif k==1
                IN(k*hexx+2*i).x=2*i;                       %8,10,12,14
                IN(k*hexx+2*i).y=k+k*sqrt(2);
                IN(k*hexx+2*i+1).x=2*i+1;                   %9,11,13,15
                IN(k*hexx+2*i+1).y=k+1+k*sqrt(2);
                IN(hexy*hexx+2*i).x=2*i+1;                %16,18,20,22
                IN(hexy*hexx+2*i).y=hexy+2*sqrt(2);
                IN(hexy*hexx+2*j+1).x=2*j+2;              %17,19,21
                IN(hexy*hexx+2*j+1).y=hexy+1+2*sqrt(2);
            else
                if mod(k,2) == 0   %number is even
                    IN(k*hexx+2*i).x=2*i;              %16,18,20,22
                    IN(k*hexx+2*i).y=k+1+k*sqrt(2);
                    IN(k*hexx+2*i+1).x=1+2*i;          %17,19,21,23
                    IN(k*hexx+2*i+1).y=k+k*sqrt(2);
                else %number is odd, so kk=3
                    IN(k*hexx+2*i).x=2*i;              %8,10,12,14
                    IN(k*hexx+2*i).y=k+k*sqrt(2);
                    IN(k*hexx+2*i+1).x=1+2*i;          %9,11,13,15
                    IN(k*hexx+2*i+1).y=k+1+k*sqrt(2);
                end
                if mod(hexy,2) == 0
                    IN(hexy*hexx+2*i).x=1+2*i;       %24,26,28,30
                    IN(hexy*hexx+2*i).y=hexy+hexy*sqrt(2);
                    IN(hexy*hexx+2*j+1).x=1+2*j+1;   %25,27,29
                    IN(hexy*hexx+2*j+1).y=hexy+1+hexy*sqrt(2);
                else
                    IN(hexy*hexx+2*i).x=2*i;         %24,26,28,30
                    IN(hexy*hexx+2*i).y=hexy+hexy*sqrt(2);
                    IN(hexy*hexx+2*j+1).x=2*j+1;     %25,27,29
                    IN(hexy*hexx+2*j+1).y=hexy+1+hexy*sqrt(2);
                end
            end
        end
    end
end

ratio=S.ElementLength/sqrt(2);

for i=1:nIN
    IN(i).x=IN(i).x*ratio;
    IN(i).y=IN(i).y*ratio;
end


%% move the x and y values in order to create a more heterogeneous network
if ~ subgeometry % this is done in main geometry otherwise
	if S.devXY>0
    	for i=1:nIN
        	IN(i).x=IN(i).x+0.3*(rand-0.5);
        	IN(i).y=IN(i).y+0.3*(rand-0.5);
    	end
	end
end
%% re-orientate the segments
% annoyingly, part of the above segments are orientated upside down. This
% gives negative flows. If the y position of second node is lower than that
% of first node, switch the nodes

if ~ subgeometry
for ie=1:nIE
	na=IE(ie).nodes(1);
	nb=IE(ie).nodes(2);
	if IN(na).y>IN(nb).y
		IE(ie).nodes=[nb na];
		% IN fields need no change, ie indicates which ie are connected but
		% not on which side
	end
end
end


%% define the lengths of the internal elements / don't change this
% calculate the lengths of the segments from the positions of the nodes,
% assuming straight segments
% might be obvious in some cases like the wheatstone, but not in others

if ~ subgeometry

	IE=LengthFromPosition(IE,IN); 
	% 2025 no longer normalized
	% % 2024 we use normalized parameters, rL as L/sum(L)
	% sumL=sum([IE.l]);
	% [IE.rL]=vout([IE.l]/sumL); 
end

%% for use as subgeometry, there are no sources or sinks, these ought to be defined in the main geometry

if subgeometry
	SE=struct; 
	nSE=0;
	[IN.nsources]=deal(0);

else
%% Define the source elements
% we will connect the high P source(s) to the bottom nodes, and the low P
% sink(s) to the top nodes in this model
% S.ConnectedBottomNodes is a vector of nodes to be connected to a high
% pressure source, counting simply from 1..ncombsx
% same for S.ConnectedTopNodes, also 1..ncombsx

% make a simple diagonal connection if not defined
if ~isfield(S,'ConnectedBottomNodes')
	S.ConnectedBottomNodes=1;
	S.ConnectedTopNodes=S.ncombsx;
end

% we need to find the indices in IN that are the bottom nodes and top nodes, 
% the bottom node indices into IN are simply 2*(1:S.ncombsx)
% the top node indices are nin+1-2*(S.ncombsx:-1:1); 

% (this was used to sort this out)
% S.ncombsx=8;S.ncombsy=5; S=honeycomb_geometry(S); yy=[S.IN.y];maxy=max(yy); find(yy>=(maxy-1e-31))

nbottom=length(S.ConnectedBottomNodes);
ntop=length(S.ConnectedTopNodes);

for se=1:nbottom
	SE(se).node=2*S.ConnectedBottomNodes(se);
end

% indices into IN of all top nodes:
v=nIN+1-2*(S.ncombsx:-1:1);
for se=1:ntop % and find the indices in v we wabt to connect
	SE(se+nbottom).node=v(S.ConnectedTopNodes(se));
end



%% define the lengths and other properties of source-connecting elements
% we don't define radii and lengths of the source elements, but rather only
% the normalized conductivity

%% define node connection table part two - the source elements
nSE=nbottom+ntop; 
% inflow and outflow connections are both 'source' but above are coupled to
% resp high and low pressure. For setting up the structure below this is
% not relevant
%no connections at start
for i=1:nIN 
    IN(i).nsources=0;
end

% fill the connections to sources
for se=1:nSE
    k=SE(se).node;
    IN(k).nsources=IN(k).nsources+1;
    IN(k).se(IN(k).nsources)=se;
end



%% define the external pressures and source/sink conductances
% note that these values might be changing in time, if they are chosen as
% input variables. In that case it should not be necessary to initialize
% them here. 

% the SE connected to bottom nodes get the high pressure, those on top the
% low pressures
[SE(1:nbottom).rPs]=deal(S.sourceP);						% (N/m2) or normalized
[SE(nbottom+1:nbottom+ntop).rPs]=deal(S.sinkP);				% (N/m2) or normalized

% and distribute the source and sink conductances equally
[SE(1:nbottom).rGs]=deal(S.sourceG/nbottom);				% [m3/s / N/m2] or normalized
[SE(nbottom+1:nbottom+ntop).rGs]=deal(S.sinkG/ntop);		% [m3/s / N/m2] or normalized

end %  if ~ subgeometry
%% current version of geometries uses IN(.).pos rather than .x and .y, such
% that we can do 3D later. Z position needs to come from main geometry, otherwise the model is 2D 
for i=1:nIN
	IN(i).pos=[IN(i).x; IN(i).y]; % have this as column vector to be able to use [IN.pos]
end
IN=rmfield(IN,{'x','y'}); % remove the x and y fields

%% assign to the output structure
S.IE=IE; S.SE=SE;S.IN=IN;
S.nIN=nIN; S.nIE=length(IE); S.nSE=nSE;


end