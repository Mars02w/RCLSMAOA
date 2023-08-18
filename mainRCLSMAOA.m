%  Hybrid Slime Mould and Arithmetic Optimization Algorithm with Random Center Learning and Restart Mutation(RCLSMAOA)
%
%  Source codes demo version 1.0                                                                      
%                                                                                                     
%  The 11th Gen Intel(R) Core(TM) i7-11700 processor with the primary frequency of 2.50GHz, 16GB memory, and the operating system of 64-bit windows 11 using matlab2021a.                                                                
%                                                                                                     
%  Author and programmer:    Hongmin Chen,Zhuo Wang,Xindong Zhou,Laith Abualigah                                                                      
%         e-Mail: jiaheminglucky99@126.com;  wang_z2002@163.com                                                                                                                                                                                                                                            


clear all 
close all
clc

N=30;  %Number of search agents
F_name='F3';     %Name of the test function
Max_iter=500;           %Maximum number of iterations

    
[lb,ub,dim,fobj]=Get_F(F_name); %Get details of the benchmark functions
[best_fun,best_position,cuve_f]=RCLSMAOA(N,Max_iter,lb,ub,dim,fobj); 


figure('Position',[454   445   694   297]);
subplot(1,2,1);
func_plot(F_name);     % Function plot
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([F_name,'( x_1 , x_2 )'])
subplot(1,2,2);       % Convergence plot
semilogy(cuve_f,'LineWidth',3)
legend('RCLSMAOA');



display(['The best-obtained solution by RCLSMAOA is : ', num2str(best_position)]);  
display(['The best optimal value of the objective funciton found by RCLSMAOA is : ', num2str(best_fun)]);  