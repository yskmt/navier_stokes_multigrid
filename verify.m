clear all
close all

P=load('points.dat');
PD=load('point_data.dat');

scatter3(P(:,1),P(:,2),P(:,3),10,PD)

PD = PD-sum(PD)/length(PD);
plot(P(:,1),PD,'-o');
hold on
plot(0:0.01:1, -sin((0:0.01:1) * 2*pi)/(2*pi)^2, '-r') 
hold off

% absolute error
sum(PD(:) + (2*pi)^-2 * sin(2*pi*P(:,1)))