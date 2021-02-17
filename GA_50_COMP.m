clear all clc clearvars
 
%%
Dataset=load('SR_div_1000.txt'); x=Dataset(:,1)';
y=Dataset(:,2)';
Min = 10;
Population = zeros(60,63); fitness=zeros(1,100000); dotplot=zeros(100000,60);

%% initial population
% generate initial population and calculate the fitness for initial fitness for i= 1:60
Population(i,:)=make_Tree(); % create 60 heap tree as initial population Y0 = tree_calculation(1,Population(i,:),x); % calculate y
Fit(i) = sum(abs(Y0(i)-y))/1000; % calculate fitness end

%% Run 100000 for each trial for z = 1:100000
dotplot(z,:)=Fit;
% to protect diversity, each new population will save 10 individual
 
 	% individuals + 50 child will become new population.
for i = 1:50
%% selection
 	% in selection function will select two pair A and B parents from
 	% top 50% of population, then two A parents and two B parents will
 	% compared with each other, next to select better one for crossover ParentA1 = Select(Fit,Population);
ParentAY1 = tree_calculation(1,ParentA1,x);

ParentA2 = Select(Fit,Population); ParentAY2 = tree_calculation(1,ParentA2,x);

%index1=randi(60); ?not used)
%ParentAY2 = tree_calculation(1,Population(index1,:),x);?not used)

%ParentA2=make_Tree();?not used)
%ParentAY2 = tree_calculation(1,ParentA2,x);?not used)

Fit_ParentA1 = sum(abs(ParentAY1-y))/1000; Fit_ParentA2 = sum(abs(ParentAY2-y))/1000;

if Fit_ParentA1 > Fit_ParentA2
%Parent A=Population (index1,:);(not used) ParentA=ParentA2;
else
ParentA=ParentA1;
end

ParentB1 = Select(Fit,Population); ParentBY1 = tree_calculation(1,ParentB1,x); ParentB2 = Select(Fit,Population); ParentBY2 = tree_calculation(1,ParentB2,x);

%index2=randi(60);?not used)
%ParentBY2 = tree_calculation(1,Population(index2,:),x);?not used)

%ParentB1=make_Tree();?not used)
%ParentBY1 = tree_calculation(1,ParentB1,x);?not used)
%ParentB2=make_Tree();?not used)
%ParentBY2 = tree_calculation(1,ParentB2,x);?not used)

Fit_ParentB1 = sum(abs(ParentBY1-y))/1000; Fit_ParentB2 = sum(abs(ParentBY2-y))/1000;

if Fit_ParentB1 > Fit_ParentB2
%Parent B=Population (index2,:);(not used) ParentB=ParentB2;
else
 
ParentB=ParentB1;
end
%% crossover
 	%create two new child
[NextGA,NextGB] = Crossover(ParentA,ParentB);
%% Mutation
 	%mutate two child
NextGA = Mutation(NextGA); NextGB = Mutation(NextGB);
 
% two parents
ParentAY = tree_calculation(1,ParentA,x); ParentBY = tree_calculation(1,ParentB,x); Fit_ParentA = sum(abs(ParentAY-y))/1000; Fit_ParentB = sum(abs(ParentBY-y))/1000;

%two childs
NextGAY = tree_calculation(1,NextGA,x); NextGBY = tree_calculation(1,NextGB,x); Fit_NextGA = sum(abs(NextGAY-y))/1000; Fit_NextGB = sum(abs(NextGBY-y))/1000;

% compare all of four individual, and find best one,
 	% if child better than parents, use child, if parent better than
 	% child, keep parent for next evolution.
if Fit_ParentA < Fit_ParentB && Fit_ParentA < Fit_NextGA && Fit_ParentA < Fit_NextGB Population(10+i,:)=ParentA;

elseif Fit_ParentB < Fit_ParentA && Fit_ParentB < Fit_NextGA && Fit_ParentB < Fit_NextGB Population(10+i,:)=ParentB;

elseif Fit_NextGA < Fit_ParentA && Fit_NextGA < Fit_ParentB && Fit_NextGA < Fit_NextGB Population(10+i,:)=NextGA;

elseif Fit_NextGB < Fit_ParentA && Fit_NextGB < Fit_ParentB && Fit_NextGA < Fit_NextGA Population(10+i,:)=NextGB;
else
Population(10+i,:)=NextGB;
end
end

%% Calculate fitness and record best one for i = 1:60
Y0 = tree_calculation(1,Population(i,:),x); Fit(i) = sum(abs(Y0-y))/1000;
end

A=sort(Fit); fitness(z) = A(1);

if fitness(z) < Min Min = fitness(z);
 



end
 
end record5(z)=Min;
DSP=[num2str(z),' : ',num2str(Min,5)]; disp(DSP);
 

%generate a tree. function Tree=make_Tree()
% heap tree structure:
% level0 1
% level1 2,3
% level2 4-7
% level3 8-15
% level4 16-31
% level5 32-63 Tree=zeros(1,63); func=[11 22 33 44 55 66];
 
% 11 means "+", 22 means"-", 33 means"x", 44 means"/", 55 means"sin" 66
% means "cos", 77 means" value of x" for i = 1:63
if i== 1 %root level R=randi(6); Tree(i)=func(R);
elseif i== 2 %1 level R2=randi(6); Tree(i)=func(R2);
elseif i== 3 %1 level R3=randi(6); Tree(i)=func(R3);
elseif rand(1)<0.5 && i>3 && i<8 %2 level R4=randi(6);
Tree(i)=func(R4);
elseif rand(1)<0.5 && i>7 && i<16 %3 level R5=randi(6);
Tree(i)=func(R5);
elseif rand(1)<0.5 && i>15 && i<32 %4 level R4=randi(6);
Tree(i)=func(R4);
else %if rand there are no node, set up as leaf nad leaf will be whether constant or x if rand(1) < 0.5
Tree(i)= 77; else
Num=linspace(-10,10,100000); index5=randi(100000); num=Num(index5); Tree(i)=num;
 

end end
 

end
 
end
 

%selection function
function Parent = Select(Fit,Population) parent=zeros(30,63);
c=sort(Fit); d=c(31:60);
get_fitness=Fit; for b=1:60
if get_fitness(b)>d(:) else
parent(b,:)=Population(b,:);
end
end
 


end
 
parent(all(parent==0,2),:)=[]; R=randi(30); Parent=parent(R,:);
 

%crossover function
function [ParentA,ParentB] = Crossover(ParentA,ParentB)
% heap tree structure:
% level0 1
% level1 2,3
% level2 4-7
% level3 8-15
% level4 16-31
% level5 32-63
% Root can not be the start point of crossover, if crossover from level 5
% it will limit the ability of crossover, which new child has less different than parents,
% diversity will lose. node in level 2 to level 4
% are best choice as start point of crossover.

%R= randi([2 31]); (note used)
%R= randi([2 15]);(note used)
%R= randi([2 8]);(note used) R= randi([4 15]);
 
%crossovver first point (start point) temp=ParentA;
ParentA(R)=ParentB(R);
ParentB(R)=temp(R);
%from start point, its subpoints will also crossover between two parents if R > 3 && R < 8
R1=2*R;
R2=2*R1;
R3=2*R2;
for i=1:2
ParentA(R1) = ParentB(R1); ParentB(R1) = temp(R1); R1=R1+1;
end
for i1 = 1:4
ParentA(R2) = ParentB(R2); ParentB(R2) = temp(R2); R2=R2+1;
end
for i2 = 1:8
ParentA(R3) = ParentB(R3); ParentB(R3) = temp(R3); R3=R3+1;
end
elseif R > 7 && R < 16 R1=2*R;
R2=2*R1;
for i=1:2
ParentA(R1) = ParentB(R1); ParentB(R1) = temp(R1); R1=R1+1;
end
for i1 = 1:4
ParentA(R2) = ParentB(R2); ParentB(R2) = temp(R2); R2=R2+1;
end end end

%mutation function
function Tree = Mutation(Tree)
% the goal of this mutation function is that to random select a point
% between level 2 to level 4, from the sub node of start point of mutation to end,
% each node or leaf will have 50% of possibility that will replace with a
% new function or constant or X.

func = [11 22 33 44 55 66 ];
% 11 means "+", 22 means"-", 33 means", 44 means"/", 55 means “sin" 66
% means "cos", 77 means" value of x"

%R = randi([2 31]);
%R = randi([8 31]);
R = randi([1 31]);

rate=0.3; % rate of mutation
% for first start point, select a new function and replace with original node. Num=randi(6);
Tree(R)=func(Num);
%Start from root if R< 2
R1=2*R;
R2=2*R1;
R3=2*R2;
R4=2*R3;
R5=2*R4;
for i=1:2
if rand(1) < rate num2=randi(6); Tree(R1) = func(num2);
 


end
 
end R1=R1+1;
 
for i1 = 1:4
if rand(1) < rate num3=randi(6);
Tree(R2) = func(num3) ;
 

end
 
end R2=R2+1;
 
for i2 = 1:8
if rand(1) < rate num4=randi(6);
Tree(R3) = func(num4) ;
 

end
 
end R3=R3+1;
 
for i3 = 1:16
if rand(1) < rate num5=randi(6); Tree(R4) = func(num5);
 

end
 
end R4=R4+1;
 
for i4 = 1:32
if rand(1) < rate Tree(R5) = 77; % x
else
r_range=10-15*rand(1); Tree(R5) = r_range;
 

end
 
end R5=R5+1;
 
elseif R > 1 && R < 4 R1=2*R;
R2=2*R1; R3=2*R2; R4=2*R3;
for i=1:2
if rand(1) < rate num2=randi(6); Tree(R1) = func(num2);
 

end
 
end R1=R1+1;
 
for i1 = 1:4
if rand(1) < rate num3=randi(6); Tree(R2) = func(num3);
 

end
 
end R2=R2+1;
 
for i2 = 1:8
if rand(1) < rate num4=randi(6); Tree(R3) = func(num4);
 

end
 
end R3=R3+1;
 
for i3 = 1:16
if rand(1) < rate Tree(R4) = 77; % x
else
r_range=10-15*rand(1); Tree(R4) = r_range;
 

end
 
end R4=R4+1;
 

elseif R >3 && R < 8 R1=2*R;
 
R2=2*R1;
R3=2*R2;
for i=1:2
if rand(1) < rate num2=randi(6); Tree(R1) = func(num2);
 

end
 
end R1=R1+1;
 
for i1 = 1:4
if rand(1) < rate num3=randi(6); Tree(R2) = func(num3);
 

end
 
end R2=R2+1;
 
for i2 = 1:8
if rand(1) < rate Tree(R3) = 77; % x
else
r_range=10-15*rand(1); Tree(R3) = r_range;
 

end
 
end R3=R3+1;
 

elseif R > 7 && R < 16 R1=2*R;
R2=2*R1;
for i=1:2
if rand(1) < rate num2=randi(6); Tree(R1) = func(num2);
end R1=R1+1;
end
for i1 = 1:4
if rand(1) < rate Tree(R2) = 77; % x
else
r_range=10-15*rand(1); Tree(R2) = r_range;
 

end else
 
end R2=R2+1;
 
R1=2*R;
for i = 1:2
if rand(1) < rate Tree(R1) = 77; % x
else
r_range=10-15*rand(1); Tree(R1) = r_range;
 


end end
 

end
 
end R1=R1+1;
 

%function for calculate the Y
%% Recursion
function Y0=tree_calculation(i,Tree,x) block=Tree(i);
if block == 11
Y0=tree_calculation(2*i,Tree,x) + tree_calculation(2*i+1,Tree,x); elseif block == 22
Y0=tree_calculation(2*i,Tree,x) - tree_calculation(2*i+1,Tree,x); elseif block == 33
 
Y0=tree_calculation(2*i,Tree,x) .* tree_calculation(2*i+1,Tree,x); elseif block == 44
Y0=tree_calculation(2*i,Tree,x) ./ tree_calculation(2*i+1,Tree,x); elseif block == 55
Y0=sin(tree_calculation(2*i,Tree,x)); elseif block == 66
Y0=cos(tree_calculation(2*i,Tree,x)); elseif block == 77
Y0= x;
else
Y0=Tree(i)*ones(1,1000);
end end
