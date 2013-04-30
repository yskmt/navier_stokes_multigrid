clear all
close all

% M=load('matrix_0.dat');
% F=load('vector_0.dat');
% scatter3(U(:,1),U(:,2),U(:,3),30,U(:,4),'fill');

M = load('test_matrix.dat');
load('test_sp_matrix.dat'); 
test_sp_matrix(:,1) = test_sp_matrix(:,1)+1;
test_sp_matrix(:,2) = test_sp_matrix(:,2)+1;
S = spconvert(test_sp_matrix);

M2 = full(S);
sum(sum(abs(M-M2)))
max_level=1;

for i=1:max_level
    U{i} = load(sprintf('results_%i.dat',i));
    legend_names{i} = sprintf('level %i', i); 
end
U{max_level+1} = load('results_100.dat');
legend_names{max_level+2} = 'final solution';
legend_names{max_level+3} = 'exact solution';

hold on
col=hsv(max_level+3);
for i=0:max_level+1
   h= plot( U{i+1}(1:end,1),U{i+1}(1:end,4), '-o');
   set(h, 'Color',col(i+1,:));
end

x=0:0.01:1;
plot(x, -(1/(2*pi))^2*sin(x*2*pi));

legend(legend_names,0);

hold off

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