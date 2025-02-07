function yMin(loc)
%
% scales x axis on current figure to place location loc at the start of the 
% figure.

curAxis = axis;
curAxis(3) = loc;
axis(curAxis);
