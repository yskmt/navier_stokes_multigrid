clear all
close all

% M=load('matrix_0.dat');
% F=load('vector_0.dat');
% scatter3(U(:,1),U(:,2),U(:,3),30,U(:,4),'fill');

max_level=2;

for i=0:max_level
    U{i+1} = load(sprintf('results_%i.dat',i));
    legend_names{i+1} = sprintf('level %i', i); 
end
U{max_level+2} = load('results_100.dat');
legend_names{max_level+2} = 'final solution';

hold on
col=hsv(max_level+2);
for i=0:max_level+1
   h= plot( U{i+1}(1:end,1),U{i+1}(1:end,4), '-o');
   set(h, 'Color',col(i+1,:));
end
hold off
legend(legend_names,0);



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