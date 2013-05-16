clear all
close all

% strong scaling
% 32^3
Ta = [0.0112691 0.00559711 0.00452089 0.00302196 0.00250888];
Tv = [0.0432389 0.041455 0.0408649 0.0634542 0.116877];
Tp = [3.77661 4.3534 4.69743 7.65514 13.1779];
Tt = [4.1107 4.66746 5.0014 7.97422 13.5518];
Nt = [16 8 4 2 1];

loglog(Nt, Ta, '-o', Nt, Tv, '-x', Nt, Tp, '-+', Nt, Tt, '-s');
xlim([0 16]);
set(gca,'XTick',[1 2 4 8 16]);
title('Strong Scaling');
legend('advection','viscosity','pressure','total');
xlabel('# of nodes');
ylabel('time [s]');

Ep = Tt(end)./(Tt.*Nt);
semilogx(Nt, Ep, '-o');
xlim([0 16]);
set(gca,'XTick',[1 2 4 8 16]);
title('Efficiency');
xlabel('# of nodes');
ylabel('efficiency');

% weak scaling
sqrt((32^3)^2/16*8); % 24*24*40
sqrt((32^3)^2/16*4); % 32*32*16
sqrt((32^3)^2/16*2); % 24*24*20
sqrt((32^3)^2/16*1); % 32*16*16

% Tw_weak = [96.2477 77.2512 70.8685 58.825 65.499];
Tw_weak = [94.9259 47.8498 19.2619 7.54821 3.63176];

semilogx(Nt, Tw_weak, '-o');
xlim([0 16]);
xlabel('# of nodes');
ylabel('total time [s]');
set(gca,'XTick',[1 2 4 8 16]);
title('Weak Scaling');
