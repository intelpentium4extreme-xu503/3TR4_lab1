function [x]=exButton( select, list)
%
% [x]=exButton( select, list)
%

if select <= length(list)
 for count = 1:length(list)
    set(list(count),'Value',0);
 end
 
 set(list(select), 'Value',1);
end
