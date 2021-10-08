% A Competitive learning-based Grey Wolf Optimizer for Engineering and its application to Power Flow Optimization Problems
%            By Aala Kalananda Vamsi Krishna Reddy and Dr.Komanapalli Venkata Lakshmi Narayana


%This is the main function
%Specify the objevtive function, population size and number of iterations here

clear all 
clc
tic
format long g;
SearchAgents_no=25; % Number of search agents

Function_name='F1'; % Name of the test function 

Max_iteration=1000; % Maximum number of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

% [Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);


[Best_score,Best_pos,GWO_cg_curve]=Clb_GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);


toc;


        



