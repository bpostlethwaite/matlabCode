%function compare(w)
wavread open8;
first_signal = ans;
vector1_length = length(first_signal);
f=0;
i=0;
s= 0;
count1 = 1;
start=0;
finish=0;
while s < 1
    if max(abs(first_signal(count1:count1+1000))) > .99
        i=i+1;
        start(i) = count1;
        f=0;
        while f < 1
            if max(abs(first_signal(count1:count1+1000))) < .99
                finish(i) = count1;
                f = 1;
            else
                count1 = count1 + 1000;
                 if count1+1000 > vector1_length
                     f = 1;
                 end
            end
        end
    end
    
    count1 = count1 + 1000;
    if count1+1000 > vector1_length
        s = 1;
    end
    
end

count_syllable = length(start);

first_syll = first_signal(start(1):finish(1));
second_syll = first_signal(start(2):finish(2));

subplot(3,1,1)
plot(first_signal)
t1=[0:1/(length(first_syll)-1):1];
subplot(3,1,2);
plot(t1,first_syll)
t2=[0:1/(length(second_syll)-1):1];
subplot(3,1,3);
plot(t2,second_syll)



length_f = length(first_syll);
length_s = length(second_syll);


fft_first = abs(fft(first_syll));
half_f = round(length_f/2);
main_f = fft_first(half_f:length_f);

fft_second = abs(fft(second_syll));
half_s = round(length_s/2);
main_s = fft_second(half_s:length_s);



t1=[0:1/(length(main_f)-1):1];
subplot(3,1,1)
plot(t1,main_f)

t2=[0:1/(length(main_s)-1):1];
subplot(3,1,2)
plot(t2,main_s)



power_f_n = sum(main_f)/length(main_f);
power_s_n = sum(main_s)/length(main_s);

avg_power_f = (abs(power_f_n-power_f)/power_f)*100;
avg_power_s = (abs(power_s_n-power_s)/power_s)*100;
if avg_power_f < 30
    if avg_power_s < 30
        
        disp('Your voice is a match')
    else
        disp('Your voice is not a match')
    end
else
    disp('Your voice is not a match')
end

            
    

