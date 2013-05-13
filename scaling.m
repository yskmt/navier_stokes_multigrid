clear all
close all

% strong scaling
% 32^3
Ta = [0.0470879, 0.0255821, 0.0162799, 0.010622 0.0148211];
Tv = [0.903047 1.22192 2.00721 3.58303 6.73668];
Tp = [191.473 279.286 462.294 846.638 1602.97];
Tt = [192.423 280.534 464.317 850.232 1609.72];
% Tw = [5.15275 6.3052 7.92871 13.0139 22.6568];
Nt = [16 8 4 2 1];
Tt = [5.05996 16; 6.12889 8; 7.5017 4; 12.67 2; 23.0299 1];

loglog(Nt, Ta, '-o', Nt, Tv, '-x', Nt, Tp, '-+', Nt, Tt, '-s');
xlim([0 16]);
set(gca,'XTick',[1 2 4 8 16]);
title('Strong Scaling');
legend('advection','viscosity','pressure','total');

Ep = Tt(end)./(Tt.*Nt);
semilogx(Nt, Ep, '-o');
xlim([0 16]);
set(gca,'XTick',[1 2 4 8 16]);
title('Efficiency');

% weak scaling
sqrt((32^3)^2/16*8); % 24*24*40
sqrt((32^3)^2/16*4); % 32*32*16
sqrt((32^3)^2/16*2); % 24*24*20
sqrt((32^3)^2/16*1); % 32*16*16

Tw_weak = [194.917 150.905 130.223 104.11 110.767];
semilogx(Nt, Tw_weak, '-o');
xlim([0 16]);
set(gca,'XTick',[1 2 4 8 16]);
title('Weak Scaling');
