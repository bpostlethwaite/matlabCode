function y = Ddiff(func,x0,h,method)
% Ddiff computes the discrete difference / derivative of a function
%
% Pass in a function func (must have one input one output,
% a value x0 to compute the difference at and
% a delta h step and a method, either 'forward' or 'center'.

fx = feval(func,x0);
fxp = feval(func,x0+h);
fxm = feval(func,x0-h);

if strcmpi(method,'forward')
    y = (fxp - fx) / h;
elseif strcmpi(method,'center')
    y = (fxp-fxm) / 2*h;
else
    disp('Please Enter either "forward" or "center"')
    y = NaN;
    return
end
    