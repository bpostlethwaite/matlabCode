% This script to put as a header in my files will load up
% often used toolboxes and function to the search path.

userdir = getenv('HOME');
f = fullfile(userdir, 'programming','matlab');

getd = @(p)path(p,path); % 

getd([f,'/toolbox_signal/']);
getd([f,'/toolbox_general/']);

