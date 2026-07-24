function S=jointrees(S1,cnv1,S2,cnv2) 
% connects two geometries, S1 and S2.
% cnv1 is a vector of nodes in S1 to be connected to the cnv2 of S2
% nv1 and nv2 must have the same size
% if there are SE connected to any nv1 or nv2, these are pruned away. 
% all fields in S1 and S2 except for the IE IN and SE and nIE nIN nSE must be
% identical
%% define the output structure
S=S1; % identical to S1; the IE IN SE nIE nIN nSE are updated below
%% cut off the sources/sinks from nodes that we are going to connect
if ~isempty(S1.SE), S1.SE=S1.SE( ~ismember( [S1.SE.node],cnv1)); end
if ~isempty(S2.SE), S2.SE=S2.SE( ~ismember( [S2.SE.node],cnv2)); end
S1.nSE=length(S1.SE);
S2.nSE=length(S2.SE);
nconnect=length(cnv1);

%% renumber the IN numbers in the IE and SE fields of S2
% renumber in the IE: nv2 will be vector for renumbering 
nv2=(1:S2.nIN)+S1.nIN; % increase the node numbers in the second tree by the number of nodes in the first tree
 for j=1:nconnect
	cn=cnv2(j); % node number to be connected in second tree
	nv2(cn)=cnv1(j);% assign node from first tree here
	nv2(cn+1:end)=nv2(cn+1:end)-1; % and lower the node numbers of the higher nodes, to prevent gaps
 end

% renumber the node numbers in the IE structure
for ie=1:S2.nIE 
	for jn=1:length(S2.IE(ie).nodes)
		S2.IE(ie).nodes(jn)=nv2(S2.IE(ie).nodes(jn));
	end
end
% renumber the node numbers in the SE structure
for se=1:S2.nSE
	S2.SE(se).node=nv2(S2.SE(se).node);
end


%% combine the IE IN and SE fields of both networks
% the fields may contain different additional information, mutually add
% missing fields and fill these with 0, assume numeric 
[S1.IE,S2.IE]=combinefields(S1.IE,S2.IE);
[S1.IN,S2.IN]=combinefields(S1.IN,S2.IN);
[S1.SE,S2.SE]=combinefields(S1.SE,S2.SE);
	function[IE1,IE2]=combinefields(IE1,IE2)
	% mutually adds fields to the both structures, such that they have the same fields. Use 0 for added values, assuming only numeric fields. 	
	fn1=fieldnames(IE1);
	fn2=fieldnames(IE2);
	IEfields=unique([fieldnames(IE1); fieldnames(IE2)]);
	newfields1=find(~ismember(IEfields,fn1));
	for k=1:length(newfields1)
		[IE1.(IEfields{newfields1(k)})]=deal(0); % generate the missing field in S1.IE and fill all with 0, assuming all fields are numeric
	end
	newfields2=find(~ismember(IEfields,fn2));
	for k=1:length(newfields2)
		[IE2.(IEfields{newfields2(k)})]=deal(0); % generate the missing field in S1.IE and fill all with 0, assuming all fields are numeric
	end
	end

S.IE=[S1.IE S2.IE];
S.nIE=S1.nIE+S2.nIE;
S.SE=[S1.SE S2.SE];
S.nSE=length(S1.SE)+length(S2.SE);

% het aantal knooppunten in de tweede staat nog op het aantal originele,
% dat zijn er nconnect te veel. We moeten uit die tweede IN nog de
% vervangen knooppunten weghalen, informatie zoals positie gaat nu wel verloren 
S2.IN=S2.IN(~ismember( 1:S2.nIN,cnv2));
S.IN=[S1.IN S2.IN]; 
% adapt the connectivity fields in IN
[S.IN,S.nIN]=MakeNodeTable(S.IE,S.SE,S.IN);
end






