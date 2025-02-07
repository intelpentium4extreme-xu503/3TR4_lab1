function [x] = mod(p,m)
%
%   mod
% 
% computes the modules of each element in p to base m
%
%
% Basic Method
%x = min(min(p));
%if (x < 0)
%  x = rem( p - floor( p/m) * m , m);
%else
%  x = rem( p , m);
%end
%
% direct method
%
x =  p - floor( p/m) * m;

