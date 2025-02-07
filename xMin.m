function xMin(loc)
%
% scales x axis on current figure to place location loc at the start of the 
% figure.

curAxis = axis;
curAxis(1) = loc;
axis(curAxis);
