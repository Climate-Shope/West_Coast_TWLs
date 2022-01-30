function chkinputdatatype(varargin)
%CHKINPUTDATATYPE Check that all inputs are double
for n = 1:nargin
    if ~isa(varargin{n},'double')
        error(message('signal:chkinputdatatype:NotSupported'));
    end
end
