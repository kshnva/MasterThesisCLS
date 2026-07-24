function IE=LengthFromPosition(IE,IN)
% calculates the length of the segments based on the position of the nodes,
% works for 2D and 3D
for i=1:length(IE)
	IE(i).l=norm(IN(IE(i).nodes(1)).pos-IN(IE(i).nodes(2)).pos);
end
