function [eq_rep_constraints,ineq_rep_constraints] = hiddenLayerConstraintsRepeated(net,u_min,u_max,repeated,u,x)

% Extract weights and biases
W = net.weights;
b = net.biases;

% Extract dimensions
dims = net.dims;
dim_hidden = dims(2:end-1);

% Extract activation function
AF = net.activation;

% Create cell for constraints
ineq_rep_constraints = {};
eq_rep_constraints = {};

% Counters
ircount = 1; %inequality
ercount = 1; %equality

% Set values for slope constraints
[Y_min,Y_max,X_min,X_max] = intervalBoundPropagation(u_min,u_max,dim_hidden,net);
if strcmp(AF, 'relu')
    alpha = 0; beta = 1;
elseif strcmp(AF, 'sigmoid')
    alpha = 0; beta = 0.25;
elseif strcmp(AF, 'tanh')
    alpha = 0; beta = 1;
end 

% All nodes connected to each other
if repeated == 1   
    % Get pre-activation variables
    for j = 1:length(dim_hidden)
        if j == 1
            x_curr_layer = x(1:dim_hidden(j));
            x_temp{j} = x_curr_layer;
            v_temp{j} = W{j}*u + b{j};
            Y_min_curr_layer{j} = Y_min(1:dim_hidden(j));
            Y_max_curr_layer{j} = Y_max(1:dim_hidden(j));
        else
            x_prev_layer = x(sum(dim_hidden(1:j-2)) + 1 : sum(dim_hidden(1:j-1)));
            x_curr_layer = x(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
            x_temp{j} = x_curr_layer;
            v_temp{j} = W{j}*x_prev_layer + b{j};
            Y_min_curr_layer{j} = Y_min(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
            Y_max_curr_layer{j} = Y_max(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
        end
    end

    for j = 1:length(v_temp)
        for jj = j:length(v_temp)
            for k = 1:length(v_temp{j})    
                if j == jj
                   for kk = k:length(v_temp{jj})
                       [alpha, beta] =  getBounds(AF,Y_min_curr_layer{j}(k),Y_max_curr_layer{j}(k),Y_min_curr_layer{jj}(kk),Y_max_curr_layer{jj}(kk));
                       alpha = 0; beta =1;
                       if j == jj && k == kk

                       else
                           x1 = x_temp{j}(k); v1 = v_temp{j}(k);
                           x2 = x_temp{jj}(kk); v2 = v_temp{jj}(kk);
                           ineq_rep_constraints{ircount} = ((x2 - x1) - alpha*(v2 - v1))*(beta*(v2 - v1) - (x2 - x1)); ircount = ircount + 1;
                       end
                   end                 
                else
                    for kk = 1:length(v_temp{jj})
                        [alpha, beta] =  getBounds(AF,Y_min_curr_layer{j}(k),Y_max_curr_layer{j}(k),Y_min_curr_layer{jj}(kk),Y_max_curr_layer{jj}(kk));
                        alpha = 0; beta =1;
                       if j == jj && k == kk

                       else
                           x1 = x_temp{j}(k); v1 = v_temp{j}(k);
                           x2 = x_temp{jj}(kk); v2 = v_temp{jj}(kk);
                           ineq_rep_constraints{ircount} = ((x2 - x1) - alpha*(v2 - v1))*(beta*(v2 - v1) - (x2 - x1)); ircount = ircount + 1;
                       end
                   end   
                end
            end
        end
    end
end

% Only nodes connected in same hidden layer
if repeated == 2
for j = 1:length(dim_hidden)
    if j == 1
        x_prev_layer = u;
        x_curr_layer = x(1:dim_hidden(j));
        v{j} = W{j}*u + b{j};
        Y_min_curr_layer = Y_min(1:dim_hidden(j));
        Y_max_curr_layer = Y_max(1:dim_hidden(j));
    else
        x_prev_layer = x(sum(dim_hidden(1:j-2)) + 1 : sum(dim_hidden(1:j-1)));
        x_curr_layer = x(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
        v{j} = W{j}*x_prev_layer + b{j};
        Y_min_curr_layer = Y_min(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
        Y_max_curr_layer = Y_max(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
    end
    for k = 1:dim_hidden(j)
        for m = k:dim_hidden(j)
            [alpha, beta] =  getBounds(AF,Y_min_curr_layer(k),Y_max_curr_layer(k),Y_min_curr_layer(m),Y_max_curr_layer(m));
            ineq_rep_constraints{ircount} = -(x_curr_layer(k) - x_curr_layer(m) - alpha)*(x_curr_layer(k) - x_curr_layer(m) - beta*(v{j}(k) - v{j}(m))); ircount = ircount + 1;  
        end
    end
end    
end

% Only nodes in adjacent hidden layers
if repeated == 3
for j = 1:length(dim_hidden)
    if j == 1
        x_prev_layer = u;
        x_curr_layer = x(1:dim_hidden(j));
        v{j} = W{j}*u + b{j};
        Y_min_curr_layer{j} = Y_min(1:dim_hidden(j));
        Y_max_curr_layer{j} = Y_max(1:dim_hidden(j));
    else
        x_prev_layer = x(sum(dim_hidden(1:j-2)) + 1 : sum(dim_hidden(1:j-1)));
        x_curr_layer = x(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
        v{j} = W{j}*x_prev_layer + b{j};
        Y_min_curr_layer{j} = Y_min(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
        Y_max_curr_layer{j} = Y_max(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
    end
    for k = 1:dim_hidden(j)
        for m = k:dim_hidden(j)
            [alpha, beta] =  getBounds(AF,Y_min_curr_layer{j}(k),Y_max_curr_layer{j}(k),Y_min_curr_layer{j}(m),Y_max_curr_layer{j}(m));
            ineq_rep_constraints{ircount} = -(x_curr_layer(k) - x_curr_layer(m) - alpha)*(x_curr_layer(k) - x_curr_layer(m) - beta*(v{j}(k) - v{j}(m))); ircount = ircount + 1;  
        end
    end
    if j > 1
        for k = 1:dim_hidden(j)
            for m = 1:dim_hidden(j-1)
                [alpha, beta] =  getBounds(AF,Y_min_curr_layer{j}(k),Y_max_curr_layer{j}(k),Y_min_curr_layer{j-1}(m),Y_max_curr_layer{j-1}(m));
                ineq_rep_constraints{ircount} = -(x_curr_layer(k) - x_prev_layer(m) - alpha)*(x_curr_layer(k) - x_prev_layer(m) - beta*(v{j}(k) - v{j-1}(m))); ircount = ircount + 1;
            end
        end
    end
end    
end

% Only for relu, enforcing equality constraints, have realised it won't do
% anything
if repeated == 4 && strcmp(AF, 'relu')
    
% Pre-processing step to find active and inactive neurons
[Y_min,Y_max,X_min,X_max] = intervalBoundPropagation(u_min,u_max,dim_hidden,net);
Ip = find(Y_min>1e-3);
In = find(Y_max<-1e-3);
Ipn = setdiff(1:sum(dim_hidden),union(Ip,In));
    
for j = 1:length(dim_hidden)
    if j == 1
        x_prev_layer = u;
        x_curr_layer = x(1:dim_hidden(j));
        v{j} = W{j}*u + b{j};
        X_min_curr_layer = X_min(1:dim_hidden(j));
        X_max_curr_layer = X_max(1:dim_hidden(j));
    else
        x_prev_layer = x(sum(dim_hidden(1:j-2)) + 1 : sum(dim_hidden(1:j-1)));
        x_curr_layer = x(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
        v{j} = W{j}*x_prev_layer + b{j};
        X_min_curr_layer = X_min(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
        X_max_curr_layer = X_max(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
    end
    for k = 1:dim_hidden(j)
        for m = k:dim_hidden(j)
            node_num_k = sum(dim_hidden(1:j-1)) + k;
            node_num_m = sum(dim_hidden(1:j-1)) + m;

            if any(node_num_k == Ip) && any(node_num_m == Ip)
                eq_rep_constraints{ercount} = x_curr_layer(k) - x_curr_layer(m) - (v{j}(k) - v{j}(m)); ercount = ercount + 1;                
            elseif any(node_num_k == In) && any(node_num_m == In)
                eq_rep_constraints{ercount} = x_curr_layer(k) - x_curr_layer(m); ercount = ercount + 1;
            else
                ineq_rep_constraints{ircount} = -(x_curr_layer(k) - x_curr_layer(m) - 0)*(x_curr_layer(k) - x_curr_layer(m) - 1*(v{j}(k) - v{j}(m))); ircount = ircount + 1;                
            end     
        end
    end
end    
end

% Only for relu, enforcing equality constraints, have realised it won't do
% anything
if repeated == 5 && strcmp(AF, 'relu')
    % Get pre-activation variables
    for j = 1:length(dim_hidden)
        if j == 1
            x_prev_layer = u;
            x_curr_layer = x(1:dim_hidden(j));
            x_temp{j} = x_curr_layer;
            v_temp{j} = W{j}*u + b{j};
        else
            x_prev_layer = x(sum(dim_hidden(1:j-2)) + 1 : sum(dim_hidden(1:j-1)));
            x_curr_layer = x(sum(dim_hidden(1:j-1)) + 1 : sum(dim_hidden(1:j)));
            x_temp{j} = x_curr_layer;
            v_temp{j} = W{j}*x_prev_layer + b{j};
        end
    end
    
    for j = 1:length(v_temp)
        for jj = j:length(v_temp)
            for k = 1:length(v_temp{j})    
                if j == jj
                   for kk = k:length(v_temp{jj})
                       if j == jj && k == kk

                       else
                           x1 = x_temp{j}(k); v1 = v_temp{j}(k);
                           x2 = x_temp{jj}(kk); v2 = v_temp{jj}(kk);
                           ineq_rep_constraints{ircount} = ((x2 - x1) - 0*(v2 - v1))*(1*(v2 - v1) - (x2 - x1)); ircount = ircount + 1;
                       end
                   end                 
                else
                    for kk = 1:length(v_temp{jj})
                       if j == jj && k == kk

                       else
                           x1 = x_temp{j}(k); v1 = v_temp{j}(k);
                           x2 = x_temp{jj}(kk); v2 = v_temp{jj}(kk);
                           ineq_rep_constraints{ircount} = ((x2 - x1) - 0*(v2 - v1))*(1*(v2 - v1) - (x2 - x1)); ircount = ircount + 1;
                       end
                   end   
                end
            end
        end
    end

end

function [alpha, beta] = getBounds(AF,Y_mini,Y_maxi,Y_minj,Y_maxj)
    Y_max2 = max(Y_maxi,Y_maxj);
    Y_min2 = min(Y_mini,Y_minj);
    if strcmp(AF, 'relu')
        if Y_max2 < 0 
            alpha = 0; beta = 0;
        elseif Y_min2 > 0 
            alpha = 1; beta = 1;
        else
            alpha = 0; beta = 1;
        end
    elseif strcmp(AF, 'sigmoid')
        %alpha = 0; beta = 0.25;
        if Y_max2 < 0 
            alpha = diffsig(Y_min2); beta = diffsig(Y_max2);
        elseif Y_min2 > 0
            alpha = diffsig(Y_max2); beta = diffsig(Y_min2);
            %alpha = min(diffsig(Y_maxi),diffsig(Y_maxj),diffsig(Y_mini),diffsig(Y_minj)); 
            %beta = max(diffsig(Y_maxi),diffsig(Y_maxj),diffsig(Y_mini),diffsig(Y_minj));
        else
            alpha = min([diffsig(Y_maxi),diffsig(Y_maxj),diffsig(Y_mini),diffsig(Y_minj)]); beta = 0.25;
        end
    elseif strcmp(AF, 'tanh')
        %alpha = 0; beta = 1;
        if Y_max2 < 0 
            alpha = difftanh(Y_min2); beta = difftanh(Y_max2);
        elseif Y_min2 > 0
            alpha = difftanh(Y_max2); beta = difftanh(Y_min2);
            %alpha = min(diffsig(Y_maxi),diffsig(Y_maxj),diffsig(Y_mini),diffsig(Y_minj)); 
            %beta = max(diffsig(Y_maxi),diffsig(Y_maxj),diffsig(Y_mini),diffsig(Y_minj));
        else
            alpha = min([difftanh(Y_maxi),difftanh(Y_maxj),difftanh(Y_mini),difftanh(Y_minj)]); beta = 1;
        end
    end    
end

function y = sig(x)
    y = 1/(1+exp(-x));
end
function y = diffsig(x)
    y = sig(x)*(1 - sig(x));
end
function y = difftanh(x)
    y = 1 - (tanh(x))^2; 
end
% if strcmp(AF, 'relu')
%         if Y_maxi < 0 && Y_maxj < 0
%             alpha = 0; beta = 0;
%         elseif Ymini >= 0 && Yminj >= 0
%             alpha = 1; beta = 1;
%         else
%             alpha = 0; beta = 1;
%         end
%     elseif strcmp(AF, 'sigmoid')
%         alpha = 0; beta = 0.25;
%         if Y_maxi < 0 && Y_maxj < 0
%             alpha = 0; beta = 0;
%         elseif Ymini >= 0 && Yminj >= 0
%             alpha = min(diffsig(Y_maxi),diffsig(Y_maxj),diffsig(Y_mini),diffsig(Y_minj)); 
%             beta = max(diffsig(Y_maxi),diffsig(Y_maxj),diffsig(Y_mini),diffsig(Y_minj));
%         else
%             alpha = min(diffsig(Y_maxi),diffsig(Y_maxj),diffsig(Y_mini),diffsig(Y_minj)); beta = 0.25;
%         end
%     elseif strcmp(AF, 'tanh')
%         alpha = 0; beta = 1;
%     end 
end