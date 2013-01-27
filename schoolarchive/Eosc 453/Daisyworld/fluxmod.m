function out = fluxmod(nflux,a,mag,itime,ntime)
%   FLUXMOD - This function returns a function which multiplies L by the specified
%   magnitude with a specified number of spikes, spaced evenly on interval.
%   nflux is number of fluxes, a is duration, mag is magnitude of
%   multiplication.

mod=ones(1,ntime);

for ii= 1:nflux

    
mod(1,ii*ntime/nflux-round(0.5*ntime/nflux)-0.5*a:ii*ntime/nflux-round(0.5*ntime/nflux))=mag;

end

out = mod(1,itime);