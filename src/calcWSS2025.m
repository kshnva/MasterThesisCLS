function S=calcWSS2025(S)
% calculates absolute WSS and delta P in S.IE from the nodal  pressures, radii and length

% somehow this code goes wrong .... 
% % % nodes=reshape([S.IE.nodes],2,[])'; % each row a pair, for each IE
% % % rP1=[S.IN(nodes(:,1)).rP];
% % % rP2=[S.IN(nodes(:,2)).rP];
% % % rL=[S.IE.rL];
% % % rWSS=S.K2*(rP1-rP2)./(2*rL);
% % % [S.IE.rWSS]=vout(rWSS);


for ie=1:S.nIE
    P1=S.IN(S.IE(ie).nodes(1)).P;
    P2=S.IN(S.IE(ie).nodes(2)).P;
    WSS_fromP(ie)=S.IE(ie).r*(P1-P2)/(2*S.IE(ie).l);
	deltaP(ie)=P1-P2;
end
[S.IE.WSS]=vout(WSS_fromP);
[S.IE.deltaP]=vout(deltaP);