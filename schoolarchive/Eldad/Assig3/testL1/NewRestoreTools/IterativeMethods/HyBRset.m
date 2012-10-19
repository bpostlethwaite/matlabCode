function options = HyBRset(varargin)
%
%   OPTIONS = HyBRset(varargin)
%
%   Create/alter options structure for HyBR code.
%   OPTIONS = HyBRset('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates an
%   options structure OPTIONS in which the named parameters have
%   the specified values.  Any unspecified parameters are set to [] (parameters
%   with value [] indicate to use the default value for that parameter when
%   OPTIONS is passed to the HyBR function). It is sufficient to type
%   only the leading characters that uniquely identify the parameter.  Case is
%   ignored for parameter names.
%   NOTE: For values that are strings, the complete string is required.
%
%   OPTIONS = HyBRset(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of OLDOPTS
%   with the named parameters altered with the specified values.
%
%   OPTIONS = HyBRset(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%   with non-empty values overwrite the corresponding old parameters in
%   OLDOPTS.
%
%   HyBRset with no input arguments and no output arguments displays all
%   parameter names and their possible values, with defaults shown in {}
%
%   OPTIONS = HyBRset(with no input arguments) creates an options structure
%   OPTIONS where all the fields are set to [].
%
%   OPTIONS = HyBRset('HyBR') creates an options structure with all
%   the parameter names and default values relevant to 'HyBR'. That is,
%           HyBRset('HyBR')
%   or
%           HyBRset(@HyBR)
%   returns an options structure containing all the parameter names and
%   default values relevant to the function 'HyBR'.
%
% HyBRset PARAMETERS for MATLAB: ( default parameter in {} )
%   InSolv - Solver for the inner problem: [ none | {Tikhonov} ]
%   RegPar - Either (a) a value for the regularization parameter
%                           [ positive scalar ]
%                   (b) a method for choosing a reg. parameter
%                           [ GCV | {modGCV} ] where
%                                   (i) 'GCV' - standard GCV
%                                   (ii) 'modGCV' - modified GCV
%                   (c) finds the optimal reg. parameter: [ optimal ]
%                                   (requires x_true)
%   Omega - If RegPar = 'modGCV', then omega is either:
%                           [ positive scalar | {adapt} ] where
%                   (a) positive scalar - a value for omega
%                   (b) adapt - uses the adaptive method
%   Iter - Maximum number of Lanczos iterations 
%                   [positive scalar | {[]} ]
%   Reorth - Reorthogonalize Lanczos subspaces: [ on | {off} ]
%   x_true - True solution : [ array | {off} ]
%                Returns error with respect to x_true at each iteration
%                and is used to compute 'optimal' regularization parameters
%   BegReg - Begin regularization after this iteration: 
%                [ positive scalar | {2} ]
%       Vx - extra space needed for finding optimal reg. parameters [{[]}]
%
%   Examples:
%     To create options with the default options for HyBR
%       options = HyBRset('HyBR');
%     To create an options structure with RegPar = 'modGCV' and Omega = .5
%       options = HyBRset('RegPar','modGCV', 'Omega',.5);
%     To change the maximum iterations to 150 in 'options'
%       options = HyBRset(options,'Iter',150);
%
%   See also HyBR.

% Print out possible values of properties. 
if (nargin == 0) && (nargout == 0)
    fprintf('            InSolv: [ none | {Tikhonov} ]\n');
    fprintf('            RegPar: [ positive scalar | GCV | {modGCV} | optimal ]\n');
    fprintf('             Omega: [ positive scalar | {adapt} ]\n');
    fprintf('              Iter: [ positive scalar  | {[]} ]\n');
    fprintf('            Reorth: [ on | {off} ]\n');
    fprintf('            x_true: [ array | {off} ]\n');
    fprintf('            BegReg: [ positive scalar | {2} ]\n');
    fprintf('                Vx: [ {[ ]} ]\n');
    return;
end

% Create a struct of all the fields with all values set to 
allfields = {'InSolv'; 'RegPar';'Omega';'Iter';'Reorth'; ...
    'x_true';'BegReg'; 'Vx'};
  
% create cell array
structinput = cell(2,length(allfields));
% fields go in first row
structinput(1,:) = allfields';
% []'s go in second row
structinput(2,:) = {[]};
% turn it into correctly ordered comma separated list and call struct
options = struct(structinput{:});

numberargs = nargin; % we might change this value, so assign it


% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
    if ischar(varargin{1})
        funcname = lower(varargin{1});
        if ~exist(funcname)
            msg = sprintf(...
                'No default options available: the function ''%s'' does not exist on the path.',funcname);
            error('Undefined function')
        end
    elseif isa(varargin{1},'function_handle')
        funcname = func2str(varargin{1});
    end
    try 
      optionsfcn = feval(varargin{1},'defaults');
    catch
        msg = sprintf(...
            'No default options available for the function ''%s''.',funcname);
        error('Undefined function', msg)
    end
    % The defaults from the HyBR functions don't include all the fields
    % typically, so run the rest of HyBRset as if called with
    % HyBRset(options,optionsfcn)
    % to get all the fields.
    varargin{1} = options;
    varargin{2} = optionsfcn;
    numberargs = 2;
end

Names = allfields;
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(['Expected argument %d to be a string parameter name ' ...
                'or an options structure \n created with HyBRset.'], i);
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                checkfield(Names{j,:},val)
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end

expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};

    if ~expectval
        if ~ischar(arg)
            error('Expected argument %d to be a string parameter name.', i);
        end

        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error( sprintf('Invalid parameter name ''%s'' ', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                error(sprintf('Ambiguous parameter name ''%s'' ', arg));
            end
        end
        expectval = 1;                      % we expect a value next

    else
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        checkfield(Names{j,:},arg);
        options.(Names{j,:}) = arg;
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error( 'Expected value for parameter ''%s''.', arg);
end


%----SUBFUNCTION---------------------------------------------
function checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'. 
%

% empty matrix is always valid
if isempty(value)
    return
end

% See if 'field' is a valid field.
validfield = true;
switch field
    case {'InSolv'} % none, tikhonov
        [validvalue, errmsg] = InSolvetype(field,value);
    case {'RegPar'} % real scalar, gcv, modgcv, optimal
        [validvalue, errmsg] = RegPartype(field,value);
    case {'Omega'} % real scalar, adapt
        [validvalue, errmsg] = Omegatype(field,value);
    case {'Iter'} % real non-negative integer
        [validvalue, errmsg] = nonNegInteger(field,value);
    case {'Reorth'} % off,on
        [validvalue, errmsg] = onOffType(field,value);
    case {'x_true'} % numeric array
        [validvalue, errmsg] = x_truetype(field,value);
    case {'BegReg'}% real non-negative integer
        [validvalue, errmsg] = nonNegInteger(field,value);
    otherwise
        validfield = false;  
        validvalue = false;
        errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
end

if validvalue 
    return;
else
  error(errmsg)
end

%------------------------------------------------------------------------

function [valid, errmsg] = InSolvetype(field,value)
% One of these strings: tikhonov
valid =  ischar(value) && any(strcmp(value,{'tikhonov','none'}));
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''tikhonov''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = RegPartype(field,value)
% One of these: positive scalar, GCV, modGCV, optimal
valid =  (isreal(value) && isscalar(value) && (value >= 0) && value == floor(value)) | (ischar(value) && any(strcmp(value,{'gcv', 'modgcv', 'optimal'})));

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a non-negative scalar or ''GCV'' or ''modGCV'' or ''optimal''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = Omegatype(field,value)
% One of these: positive scalar, numeric array, adapt
valid =  (isreal(value) && isscalar(value) && (value >= 0) )|(isreal(value) && isnumeric(value) && all(value>0) ) | (ischar(value) && any(strcmp(value,{'adapt'})));

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a non-negative scalar or ''adapt''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = nonNegInteger(field,value)
% Any nonnegative real integer scalar
valid =  isreal(value) && isscalar(value) && (value >= 0) && value == floor(value) ;

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative integer.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = onOffType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'on';'off'}));
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = x_truetype(field,value)
% Either a numeric array or off
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'off'})));
if ~valid
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array.',field);
else
  errmsg = '';
end
