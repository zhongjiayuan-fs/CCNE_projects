
clear;
clc;
load pericyte_to_neuron.mat;
for t = 1:6
    time_entropy = entropy_matrix(:,:,t);
    for c = 1:cell_num(t)
        [sorted_H,idx] = sort(time_entropy(:,c),'descend');
        SH(t,c) = sum(sorted_H(1:total_node_num*0.05));
    end
    result(t) = mean(SH(t,1:cell_num(t)));
end
figure;
combineData = [SH(1,1:cell_num(1)),SH(2,1:cell_num(2)),SH(3,1:cell_num(3)),SH(4,1:cell_num(4)),SH(5,1:cell_num(5)),SH(6,1:cell_num(6))];  
group = [2*ones(1,cell_num(1)),3*ones(1,cell_num(2)),4*ones(1,cell_num(3)),5*ones(1,cell_num(4)),6*ones(1,cell_num(5)),7*ones(1,cell_num(6))]; 
boxplot(combineData, group, 'notch', 'on', 'sym', 'r+', 'Colors', [0.9, 0.1, 0.1]);
hold on;
t=1:6;
plot(t, result, 'Color', [0.9, 0.1, 0.1], 'LineWidth', 3);
hold on;
scatter(1:6, result, 100, [0.9, 0.1, 0.1], 'filled');
ylabel( 'CSCNE');
ylim([0,40]); 
yticks(0:10:40); 
box off








