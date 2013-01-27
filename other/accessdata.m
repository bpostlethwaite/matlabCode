% Access Database

addpath sac
addpath Data
addpath Functions
user = getenv('USER');
datadir = ['/home/',user,'/Dropbox/ComLinks/Programming/matlab/thesis/Data'];
databasedir = [datadir,'/database'];
load(sprintf('%s/database.mat',datadir))


