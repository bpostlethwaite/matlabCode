clear all; close all;
n = 256;

cmaps = {'jet', 'hsv' ,'hot', 'cool', 'spring', 'summer', 'autumn',...
    'winter', 'gray', 'bone', 'copper'};  %,'pink', 'lines'};


div = 1:n;
ax = linspace(0, 1, n);






for ii = 1:length(cmaps)
    m(ii).type = cmaps{ii};
    fh = str2func(m(ii).type);
    m(ii).map = fh(n);
    
    I = floor(abs(diff(diff(m(ii).map))) * 1000) > 0;

    Ir = [1, div(I(:,1)) + 2, n];
    Ig = [1, div(I(:,2)) + 2, n];
    Ib = [1, div(I(:,3)) + 2, n];
    
    %figure()
    %plotmap(m(ii), Ir, Ig, Ib)
    
    
    rdiv = [sprintf('%1.3f, ', ax(Ir))];
    rdiv = ['[', rdiv(1:end - 2), ']'];
    
    rval = [sprintf('%1.3f, ', m(ii).map( Ir, 1)) ];
    rval = ['[', rval(1:end - 2), ']'];
    
    r = sprintf('  r: [\n    %s\n  , %s\n ]\n', rdiv, rval);
    
    gdiv = [sprintf('%1.3f, ', ax(Ig))];
    gdiv = ['[', gdiv(1:end - 2), ']'];
    
    gval = [sprintf('%1.3f, ', m(ii).map( Ig, 2)) ];
    gval = ['[', gval(1:end - 2), ']'];
    
    g = sprintf(', g: [\n    %s\n  , %s\n ]\n', gdiv, gval);
    
    bdiv = [sprintf('%1.3f, ', ax(Ib))];
    bdiv = ['[', bdiv(1:end - 2), ']'];
    
    bval = [sprintf('%1.3f, ', m(ii).map( Ib, 3)) ];
    bval = ['[', bval(1:end - 2), ']'];
    
    b = sprintf(', b: [\n    %s\n  , %s\n ]\n', bdiv, bval);
    
    
    
    fprintf(', %s: {\n%s%s%s\n  }\n', m(ii).type, r,g,b)
    
    
end