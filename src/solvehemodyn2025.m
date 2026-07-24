function [S]=solvehemodyn2025(S,modelpart)
% usage: [S]=solvehemodyn2025(S);
% version 2025 using ABSOLUTE parameters
% from version 20 october 2016 Janina Schwarz
% This solves the network pressures and flows, requiring network
% topological info and SE(.).Gs Ps , IE(.).G
% it returns SE(.).Qs , IE(.).Q IE(.).P (midway pressure), IN(.).P
% (nodal pressure)

% 2024 added modelpart ('electrical' or 'hemodynamic') to work on either
% system part

%% Usage
% This function solves the flow and pressure distribution throughout a network
% defined by nodes and elements.
%
% * IN: array of struct for internal nodes, requires fields _cn_ (connected nodes),
% _ie_ and _se_
% * IE: array of stuct for internal elements, requires fields  _rG_ (relative conductance) and
% _nodes_
% * SE: array of struct for source elements, requires fields _rGs_ (conductance), _rPs_ and
% _node_
%
% S.I. units are used: rradius (_r_) in [-], rpressure (_P_) in [-], relative flow
% (_rQ_) in [-], relative conductance (_rG_) in [-] and relative wall shear
% stress (rWSS) in [-]

%% Background
% At every node $j$, net inflow equals net outflow, i.e.
%
% $$\sum_{n\neq j} (P_j-P_n)*G_{jn} = 0$$
% 
% where $P_n$ is the pressure at node $n$ and
% $G_{jn}$ is the conductance of the element connecting nodes $j$ and $n$.
% All intenal nodes together yield a system of linear equations. This
% linear system can be rewritten in matrix notation with split
% contributions from connected internal nodes (at unknown pressure _P_) 
% and source elements (with source nodes at known pressure _Ps_):
% _A*P+B*Ps=0_.
%% define whether this is hemodynamic or electrical part
% for solving the electrical  system we temporarily paste that into the
% hemodynamic fields

if nargin==1, modelpart='hemodynamic'; end

if strcmp(modelpart,'electrical')
    % switch fields and switch them back at end
    holdIE=S.IE; holdnIE=S.nIE; S.IE=S.EIE; S.nIE=S.nEIE;
    holdIN=S.IN; holdnIN=S.nIN; S.IN=S.EIN; S.nIN=S.nEIN;
    holdSE=S.SE; holdnSE=S.nSE; S.SE=S.ESE; S.nSE=S.nESE;
end


%% Definition of A, B & Ps
% A is a system matrix defining internal connectivity.
% Diagonals reflect currents flowing away from the node, other elements currents
% flowing towards nodes.
%
% B is a matrix defining connectivity to sources. B has the same number of
% rows as A and one column per pressure source, each element is a
% conductance to a pressure source.
% [js]: Please note that I'm using a slightly different definition of A & B
% than in _solvehemodyn.m_.
% This definition also works for multiple source elements connected to one
% node and for multiple elements connecting the same two nodes.

nIN=length(S.IN); nIE=length(S.IE); nSE=length(S.SE);
if ~isfield(S.IN(1),'nconnect')
    for i=1:nIN, S.IN(i).nconnect=length(S.IN(i).cn); S.IN(i).nsources=length(S.IN(i).se); end
end

% find parallel IEs
IEnodes=sort(reshape([S.IE.nodes],[2,nIE])',2);

% define A & B as sparse matrices
cA=1;cB=1;

% all INs
for j=1:nIN
    cAs=cA;cBs=cB; % first indices for this row
    cn=unique(S.IN(j).cn);
    cIEs=S.IN(j).ie;
    IEnodesthisIN=IEnodes(cIEs,:);
    % do something similar also for SEs?
    for i=1:length(cn)
        conIE=cIEs(logical(sum(ismember(IEnodesthisIN,cn(i)),2)));
        % computation of joint parallel conductance
        G=0;
        for k=1:length(conIE), G=G+S.IE(conIE(k)).G; end
        vA(cA)=G;    
        colA(cA)=cn(i);
        cA=cA+1;
    end
    d=-sum(vA(cAs:cA-1));
    % then define B
    for i=1:S.IN(j).nsources
        vB(cB)=S.SE(S.IN(j).se(i)).Gs;
        colB(cB)=S.IN(j).se(i);
        d=d-vB(cB);
        cB=cB+1;
    end
    % diagonal term of A
    vA(cA)=d;
    colA(cA)=j;
    rowA(cAs:cA)=j;rowB(cBs:cB-1)=j; % one row for every IN
    cA=cA+1;
end

A=sparse(rowA,colA,vA,nIN,nIN);
B=sparse(rowB,colB,vB,nIN,nSE);
clear v* col* row* IEnodes % save some memory 

% vector with source pressures
PS=sparse([S.SE.Ps]');

% full(A)
% full(B)
% full (PS)
% BPS=(B*PS); full(BPS)


%% Solution by matrix inversion
% Find pressures at internal nodes and in elements:
 
% spparms('spumoni',1); % info about solver
P=-A\(B*PS); % actual solution by matrix inversion
P=full(P);

% distribute over the elements
for j=1:nIN, S.IN(j).P=P(j); end  

% flows and mean pressures in elements
for i=1:nIE
    P1 = S.IN(S.IE(i).nodes(1)).P;      %Pressure node 1 van het ith interne element
    P2 = S.IN(S.IE(i).nodes(2)).P;      %Pressure node 2 van het ith interne element
    S.IE(i).Q = (P1-P2) * S.IE(i).G;   
    S.IE(i).P = 0.5 * (P1+P2);    
end

% for source elements, flow is positive into the network, venous flow thus
% shows as negative
for i=1:nSE
    P1 = S.SE(i).Ps; 
    P2 = S.IN(S.SE(i).node).P;
    S.SE(i).Qs = (P1-P2) * S.SE(i).Gs; 
    %SE(i).P = 0.5 * (P1+P2);
end

%% define whether this is hemodynamic or electrical part
% for solving the electrical  system we temporarily paste that into the
% hemodynamic fields

if strcmp(modelpart,'electrical')
    % switch fields and switch them back at end
    S.EIE=S.IE; 
    S.nEIE=S.nIE;
    S.IE=holdIE; 
    S.nIE=holdnIE; 
    
    S.EIN=S.IN; 
    S.nEIN=S.nIN;
    S.IN= holdIN; 
    S.nIN=holdnIN; 
    
    S.ESE= S.SE; 
    S.nESE=S.nSE;
    S.SE=holdSE; 
    S.nSE=holdnSE; 
    
end