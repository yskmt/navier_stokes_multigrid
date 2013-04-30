clear all
close all


T = [5.2855 7.68844 11.5905 21.6132 41.5713];
P = [16 8 4 2 1];

loglog(P,T,'-o');
xlabel('# of threads');
ylabel('time [s]');
title('64^3 degree of freedoms')
xlim([1 16]);

(log(T(1))-log(T(5))) / (log(P(1))-log(P(5)))

T(5)/T(1)

T = [5.73453  2.58304 1.26081 0.701893 0.50098]
P = [16 8 4 2 1];
loglog(P,1./sqrt(T),'-o');