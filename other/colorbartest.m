% colorbar test
close all;
clear all;

J = jet(81);

for ii = 1:81
    fprintf(', "#%s"\n', rgbconv(J(ii,:)) );
end