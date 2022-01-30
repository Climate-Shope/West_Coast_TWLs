function [h, errid, errmsg, msgObj] = uddpvparse(c, varargin)
%UDDPVPARSE   Parse p/v pairs using a udd object.
%   UDDPVPARSE(C) Returns a UDD object of class c.
%
%   UDDPVPARSE(C, H) Returns H if it is of class C.  Errors if it is not
%   the correct class.
%
%   UDDPVPARSE(C, P1, V1, P2, V2, etc.) Returns a UDD object of class C
%   with its parameters P1, P2, etc. set to V1, V2, etc.
%
%   UDDPVPARSE(C, CELLC, ...) If the second argument is a cell array, it
%   will be used to construct the default object.
%
%   [H, ERRID, ERRMSG] = UDDPVPARSE(...) Returns an error identifier in
%   ERRID and an error message in ERRMSG. These will be empty when there
%   are no errors.  The returned ERRID and ERRMSG can later be passed
%   directly to ERROR.  This is used to avoid deep stack traces.  If one H
%   is requested, UDDPVPARSE will error when ERRMSG is not empty.
%
%   [H, ERRID, ERRMSG, MSGOBJ] = UDDPVPARSE(...) Returns a message object
%   that can be used to call error(MSGOBJ). MSGOBJ is empty when no error
%   occurred.
%
%   This function is meant to be called by functions which can take either
%   an options object (C) or P/V pairs for the options object.
%
%   Examples:
%
%   % Constructs a default object
%   uddpvparse('fdopts.sosscaling') 
%
%   % Uses the input as the object
%   h = fdopts.sosscaling;
%   h.MaxNumerator = 10;
%   uddpvparse('fdopts.sosscaling', h)
%
%   % Uses the P/V pairs
%   uddpvparse('fdopts.sosscaling', 'MaxNumerator', 10)
%
%   % Uses the optional constructor syntax
%   Hd = dfilt.df2sos;
%   set(Hd, 'Arithmetic', 'Fixed', ...
%       'CoeffAutoScale', false, ...
%       'NumFracLength', 17);
%   uddpvparse('fdopts.sosscaling', {'scaleopts', Hd}, ...
%       'ScaleValueConstraint', 'po2')
error(nargchk(1,inf,nargin,'struct'));

identifier = '';

h = [];

if nargin == 1,
    h = feval(c);
else

    % If the first extra is a cell array, it must be the "constructor".
    if iscell(varargin{1})
        cons = varargin{1};
        varargin(1) = [];
    else
        cons = {c};
    end

    if length(varargin) == 0
        try
            h = feval(cons{:});
        catch ME
            throw(ME);            
        end
    else
        % Check if the first input is actually the object.  Cannot use ISA
        % because that returns true for subclasses.
        if strcmpi(class(varargin{1}), c)

            % If the object is the input, just return it and ignore the rest of
            % the inputs.
            h = varargin{1};
            varargin(1) = [];
        end
        
        if ~isempty(varargin)
          
            if rem(length(varargin), 2)

                % Cannot have an odd # of parameters
                identifier = 'signal:uddpvparse:invalidPVPairs';
            else

                % Assume that the rest of the inputs must be P/V pairs.  If they
                % are not they will be captured in the 'try/catch'.
                if isempty(h)
                    h = feval(cons{:});
                end
                try
                    set(h, varargin{:});
                catch ME
                    throw(ME);
                end
            end
        end
    end
end

if ~isempty(identifier)
    errid = identifier;
    msgObj = message(identifier);
    errmsg = getString(msgObj);    
else
    errid = '';
    errmsg = '';
    msgObj = [];
end

% Only error if the msg structure is not requested.
if nargout == 1
    if ~isempty(errmsg), error(msgObj); end
end