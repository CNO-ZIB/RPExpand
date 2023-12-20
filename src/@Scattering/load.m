function varargout = load(filename)
%LOAD load cartesian export with derivatives for the point evaluation
%
%   [out,names] = LOAD(sc,filename) returns an (n,3) matrix where n is the
%   number of partial derivatives and the last axis accounts for x, y and z 
%   direction. The names refer to the derivative parameters.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

id = fopen(filename,'r');
ndx = 1; name = cell(1,10);
while 1
    tline = fgetl(id);
    if ~ischar(tline), break, end
    if startsWith(tline,'ParameterDerivatives')
        ndx = ndx+1;
        name{ndx} = ['d_' regexp(tline,'(?<==)\w*','match','once')];
    elseif startsWith(tline,'TensorQuantityVector:1:Type')
        name{1} = regexp(tline,'(?<==)\w*','match','once'); ndx = 1;
    elseif startsWith(tline,'#Field')
        s.(name{ndx}) = zeros(3,1); 
        for it = 1:3
            rpart = regexp(fgetl(id),'-?\d+\.\d+e[+-]\d+','match','once');
            ipart = regexp(fgetl(id),'-?\d+\.\d+e[+-]\d+','match','once');
            s.(name{ndx})(it) = str2double(rpart)+1i*str2double(ipart);
        end
        ndx = ndx+1;
    end
end
fclose(id);
varargout{1} = struct2array(s).';
if nargout>1, varargout{2} = fieldnames(s); end
end
