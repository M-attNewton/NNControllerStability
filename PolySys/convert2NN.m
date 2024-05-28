load('nn_3_relu_data.mat')
data = nn_3_relu_data;
dim_in = data(1);
dim_out = data(2);

NN_data = data(6:end);

W{1} = zeros(20,dim_in); 
b{1} = zeros(20,1);

W{2} = zeros(20,20); 
b{2} = zeros(20,1);

W{3} = zeros(1,20); 
b{3} = zeros(1,1);

k = 1;
for ii = 1:20
    for jj = 1:dim_in
        W{1}(ii,jj) = NN_data(k);
        k = k + 1;
    end
    b{1}(ii) = NN_data(k);
    k = k + 1;
end

for ii = 1:20
    for jj = 1:20
        W{2}(ii,jj) = NN_data(k);
        k = k + 1;
    end
    b{2}(ii) = NN_data(k);
    k = k + 1;
end

for ii = 1:1
    for jj = 1:20
        W{3}(ii,jj) = NN_data(k);
        k = k + 1;
    end
    b{3}(ii) = NN_data(k);
    k = k + 1;
end

return

%% another function for julia script

dims = [2,3*ones(1,100),1];
start_nodes = 6;
num_nodes = 8;
W{1} = W{1}(start_nodes:num_nodes,:);
W{2} = W{2}(start_nodes:num_nodes,start_nodes:num_nodes);
W{3} = W{3}(1,start_nodes:num_nodes);
b{1} = b{1}(start_nodes:num_nodes);
b{2} = b{2}(start_nodes:num_nodes);
b{3} = b{3}(1);

AF = 'relu';
net.activation = 'relu';
net.dims = dims;
net.weights{1} = W{1};
net.biases{1} = b{1};

net.weights{101} = W{3}
net.biases{101} = b{3};

for j = 2:100
    net.weights{j} = W{2};
    net.biases{j} = b{2};
end

for i = 1:(length(dims) - 1)
    weights(1:dims(i+1), 1:dims(i), i) = net.weights{i};
    biases(1:dims(i+1), i) = net.biases{i};
end

return
filename = 'ReachSparsePsatz/newSys2D_3by100';
save(filename,'weights','biases','dims','AF')
save(filenameNET,'net')