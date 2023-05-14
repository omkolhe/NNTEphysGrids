function [vx, vy] = wavefront_direction(pd,s)
% *WAVE*
%
% WAVEFRONT DIRECTION    calculates the vecocity vector 
%
% INPUT
% pd - phase gradient direction
% s -  instantaneous speed 
%
% OUTPUT
% vx - the x component of the vector pointing in velocity direction
% vy - the y component of the vector pointing in velocity direction
%

assert( ( ndims(pd) == 3 ) & ( ndims(s) == 3 ), 'datacube inputs required' );
assert( isequal( size(pd), size(s) ), 'datacube sizes must be equal' );

mag = -1;
% mag = -abs(s)/10; %in cm/s
vx = squeeze(mean((cos(pd).*mag),[1 2]))' ;
vy = squeeze(mean((sin(pd).*mag),[1 2]))';
end

