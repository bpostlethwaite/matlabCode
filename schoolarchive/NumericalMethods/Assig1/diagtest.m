

s = 1:10;
s = s'

c = ones(length(s),1);

spp = spdiags([s,-s],[-1,1],length(s),length(s));
spy(spp)
spp = spp'
full(spp)
