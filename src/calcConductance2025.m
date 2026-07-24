function [S]=calcConductance2025(S)
% calculates absolute (not relative as in loop paper) conductance in the IE of S
G=pi*[S.IE.r].^4./(8*S.viscosity*[S.IE.l]);
[S.IE.G]=vout(G);