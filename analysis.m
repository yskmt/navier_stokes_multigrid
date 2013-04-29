clear all
close all

M=load('matrix_3.dat');
F=load('vector_3.dat');
U = load('results.dat');

% spy(matrix);

[n,m]=size(M);


% 
% scatter3(U(:,1),U(:,2),U(:,3),30,U(:,4),'fill');
% plot(U(:,1),U(:,4), '-o');
% hold on
% plot(0:0.01:1, -sin((0:0.01:1) * 2*pi)/(2*pi)^2, '-r') 
% hold off


% v=zeros(9,9,9);

% for i=1:length(U)
%     v(U(i,1)+1,U(i,2)+1,U(i,3)+1) = U(i,4);
% end
%     
% x=U(:,1);
% y=U(:,2);
% z=U(:,3);
% % v=U(:,4);
% 
% [x y z v] = flow;
% h=contourslice(x,y,z,v,[1:9],[],[0], linspace(-8,2,10));
% axis([0 10 -3 3 -3 3]); daspect([1 1 1])
% camva(24); camproj perspective;
% campos([-3 -15 5])
% set(gcf, 'Color', [.3 .3 .3], 'renderer', 'zbuffer')
% set(gca, 'Color', 'black' , 'XColor', 'white', ...
%                'YColor', 'white' , 'ZColor', 'white')
% box on