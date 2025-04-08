
clear;
clc;
close all;
[mydata,pipi]= xlsread('pericyte_to_neuron_data.xlsx');
mydata=mydata(:,2:624);


mydata=log(mydata+1);
fid=fopen('Undirected_ppi_network.txt');
adjacent_network={};
j=0;
while ~feof(fid)
    tline=fgetl(fid);
    j=j+1;
    adjacent_network{j}=regexp(tline, '\t', 'split');
end
fclose(fid);
total_node_num=j;
cell_num=[76,86,48,283,61,69]; 
time_point=6;


time_1_data=mydata(:,1:cell_num(1));
time_2_data=mydata(:,1+cell_num(1):sum(cell_num(1:2)));
time_3_data=mydata(:,1+sum(cell_num(1:2)):sum(cell_num(1:3)));
time_4_data=mydata(:,1+sum(cell_num(1:3)):sum(cell_num(1:4)));
time_5_data=mydata(:,1+sum(cell_num(1:4)):sum(cell_num(1:5)));
time_6_data=mydata(:,1+sum(cell_num(1:5)):sum(cell_num(1:6)));

nn=1;
p_value=0.05;
count=2;
mix_numcell=min(cell_num);
eps=1e-10;
nn=1;
p_value=0.05;
for t=1:time_point
    input_data=mydata(:,nn:sum(cell_num(1:t)));
    nn=1+sum(cell_num(1:t));
    input_data_size=size(input_data);
    k=round(0.1*input_data_size(2));      
    for na=1:total_node_num
        center=adjacent_network{na}{1};
        num=0;
        clear localnet_x2yMI_weigtht_matrix;
        clear localnet_x2y_cond_matrix;
        clear localnet_y2xMI_weigtht_matrix;
        clear localnet_y2x_cond_matrix;
        clear localnet_weigtht_matrix;
        for n=2:length(adjacent_network{na})
            nei=adjacent_network{na}{n};
            x = input_data(str2num(center),:);
            y = input_data(str2num(nei),:);
            out_x2y = MI_cnet(x',y',k);
            num=num+1;
            localnet_weigtht_matrix(num,1:cell_num(t))=out_x2y;
        end
       normalized_matrix = localnet_weigtht_matrix ./ sum(localnet_weigtht_matrix);
       normalized_matrix = fillmissing(normalized_matrix,'constant',0);
        
        num_cols = size(localnet_weigtht_matrix, 2); %            
        column_entropies = zeros(1, num_cols);
        for col = 1:num_cols
            out_column_data = normalized_matrix(:, col);
            outnum=nnz(normalized_matrix(:, col));
            
            center_data=input_data(str2num(center),:);
            sd_center=std(center_data);
            center_data(col)=[];
            sd_center_remove_col=std(center_data);
            sd_delt=sqrt(input_data_size(2))*abs(sd_center-sd_center_remove_col);
            H = -sum(out_column_data .* log2(out_column_data + eps));
            column_entropies(col) =(1/(outnum))*H*sd_delt;
        end
        column_entropies = fillmissing(column_entropies,'constant',0);
        entropy_matrix(na,1:num_cols,t)=column_entropies;
    end
    t
end
entropy_matrix = fillmissing(entropy_matrix,'constant',0);
save pericyte_to_neuron.mat;








