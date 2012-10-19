% Sound Gen

dw   = 0.02;
w    = 066.16;
w1   = w - dw*w;
w2   = w + dw*w;
tmax = 1;
t = [0:0.000125:tmax];
v = sin(2*pi*w1*t) + sin(2*pi*w2*t);

wavwrite(v, 'asharp.wav');


f = sin(2*pi*174.61*t);
g = sin(2*pi*195.99*t);
a = sin(2*pi*220*t);
b = sin(2*pi*123.47*t);
c = sin(2*pi*130.81*t);
d = sin(2*pi*146.83*t);
e = sin(2*pi*164.81*t);
 
%Now to create a line of music use the following command:
 
line1 = [c,c,b,c,b,c];
line2 = [e,f,d,f,d,e];
line3 = [g,a,g,a,g,g];
line4 = [v,v,v,v,v,v];
%The letters should represent the notes that you have created in MatLab.  Put the notes in the order you want them to play.
%To create your song use:
 
song = [line1 + line2+ line3 + 0.4*line4];

wavwrite(song, 'mysong.wav')