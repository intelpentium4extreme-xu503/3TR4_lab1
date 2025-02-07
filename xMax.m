function xMax(loc)
%
% scales x axis on current figure to place location loc at the start of the 
% figure.

curAxis = axis;
curAxis(2) = loc;
axis(curAxis);
