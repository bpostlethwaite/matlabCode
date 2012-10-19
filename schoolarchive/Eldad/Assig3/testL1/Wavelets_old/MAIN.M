% the main function of the project! 
% note: auread the sound file first

function y = main(x)

x = normalize(x);
x=x';

x_padded = zero_pad(x);

le=32;
h=daubcqf(le,'min');
L=7;

[x_padded, L] = mdwt(x_padded, h, L);

%load OTTFF

L3_x = x_padded(1025:2048);  	% consider the 3rd octave
U = mean(abs(L3_x))/5;

if(U >= .0010)  	% the number is 1, 4, or 5

	L4_x = x_padded(513:813); 	% consider the first part of the 4th octave
	s = find(abs(L4_x/5.3) > .005); 
	S=length(s);

	if(S >= 100)
		y = 4;
	
	elseif(S <= 80 & S >= 50)
		y = 1;

	elseif(S <= 35)
		y = 5;

	else	
	one_diff = normdiff(ONE,x_padded);
	four_diff = normdiff(FOUR,x_padded);
	five_diff = normdiff(FIVE,x_padded);
	diff1 = [one_diff four_diff five_diff];
	M = max(diff1);

		if (M == one_diff)
			y = 1;
		elseif (M == four_diff);
			y = 4;
		else
			y = 5;
		end
	end


else	 	% the number is 2 or 3

	L2 = x_padded(2049:4096);	% look at the 2nd octave 
	L2 = L2/abs(max(L2)); 
	a = find(abs(L2) > .4);
	A=length(a);

	if(A >= 90)
		y = 2;

	elseif(A < 50)
		y = 3;

	else
	two_diff = normdiff(TWO,x_padded);
	three_diff = normdiff(THREE,x_padded);
	diff2 = [two_diff three_diff];
	M = max(diff2);

		if (M == two_diff)
			y = 2;
		else 
			y = 3;
		end
	
	end
end


fprintf('\n The number is   %d! \n\n', y)












