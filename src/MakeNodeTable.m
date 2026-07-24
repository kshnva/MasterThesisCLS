function [IN,nin]=MakeNodeTable(IE,SE,IN)
% 18 feb 2022: added IN as optional input, this maintains existing fields,
% not related to connectivity
%% DEFINITION OF THE CONNECTIVITY OF THE MODEL: the node table
% now generate a table of internal nodes from the info you have above
% this works for nodes with as many connections as you want, but typically
% 3
%==
% added october 2022 for FP2P backwards compatible with old version
if nargin<3, IN=[]; end
%==
nie=length(IE);
nse=length(SE);

% total number of internal nodes
nin=max([IE.nodes]);  

%no connections at start
for i=1:nin 
    IN(i).nconnect=0;
    IN(i).nsources=0;
	% added 23/3/22 to make sure old connectivity data disappear
	IN(i).cn=[];
	IN(i).ie=[];
	IN(i).se=[];
end

% fill the internal connections based on the definition of IE,  
for ie=1:nie
    node(1)=IE(ie).nodes(1); 
    node(2)=IE(ie).nodes(2);
    %G=IE(ie).G; kan hier weg?
    for j=1:2
        k=node(j);
        IN(k).nconnect=IN(k).nconnect+1;
        nc=IN(k).nconnect;
        IN(k).cn(nc)=node(3-j); % 3-j means 1=>2 and 2=>1
        IN(k).ie(nc)=ie; % points back to element table
    end
end

% fill the connections to sources
for se=1:nse
    k=SE(se).node;
    IN(k).nsources=IN(k).nsources+1;
    IN(k).se(IN(k).nsources)=se;
end
