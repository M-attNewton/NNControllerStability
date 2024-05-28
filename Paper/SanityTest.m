%% test script

syms t1 t2

load('inv_pend_MPC_bias_free.mat')
mass = 0.15;
leng = 0.5;
mu = 0.5;
grav = 9.81;

%u = W{3}*(tanh(W{2}*(tanh(W{1}*[t1;t2] + b{1})) + b{2})) + b{3};
%u = W{3}*(tanh(W{2}*(tanh(W{1}*[t1;t2]))));
u = max(-1,min(W{3}*(tanh(W{2}*(tanh(W{1}*[t1;t2])))), 1));
testsol = 429.1822*t1^4 + 42.042*t1^3*t2 + 24.724*t1^2*t2^2 + 3.3806*t1*t2^3 + 0.39297*t2^4 + 9.7491e-12*t1^3 + 4.856e-13*t1^2*t2 + 2.0667e-13*t1*t2^2 + 1.1171e-14*t2^3 + 0.0021561*t1^2 + 0.00015804*t1*t2 + 4.365e-05*t2^2
 %u = 1
testdot1 = t2;
testdot2 = (mass*grav*leng*sin(t1) - mu*t2 + u)/(mass*leng^2);

testderV = diff(testsol,t1)*testdot1 + diff(testsol,t2)*testdot2;
 
fsurf(testderV,[-0.3,0.3,-1.4,1.4]) %[-10,10,-10,10])%,[-0.3,0.3,-1.4,1.4])
hold on
fsurf(0)