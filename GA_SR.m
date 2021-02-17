clear all clc clearvars
%% initialization Dataset=load('SR_div_1000.txt'); x=Dataset(:,1)';
y=Dataset(:,2)';
%% set initial Min fitness as 10 Min = 10;
Population = zeros(60,63); fitness=zeros(1,100000); dotplot=zeros(100000,60); P1=zeros(1,1000); P2=zeros(1,1000);
%% initial population
% generate initial population and calculate the fitness for initial fitness for i= 1:60
Population(i,:)=make_Tree(); % create 60 heap tree as initial population Y0 = tree_calculation(1,Population(i,:),x); % calculate y
Fit(i) = sum(abs(Y0(i)-y))/1000; % calculate fitness end

%% Run 100000 for each trial for z = 1:100000
dotplot(z,:)=Fit;
% to protect diversity, each new population will save 10 individual

 	% individuals + 50 childs will become new population. for i = 1:50
%% selection
 	% in selection function will select two pair A and B parents from
 	% top 50% of population, then two A parents and two B parents will
 	% compared with each other, next to select better one for crossover A=randi(60);
ParentA1=Population(A,:); P1=tree_calculation(1,ParentA1,x);

B=randi(60); ParentA2=Population(B,:); P2=tree_calculation(1,ParentA2,x);

Fit_ParentA1 = sum(abs(P1-y))/1000; Fit_ParentA2 = sum(abs(P2-y))/1000;

if Fit_ParentA1 > Fit_ParentA2
%ParentA=Population(index1,:); ParentA=ParentA2;
else
ParentA=ParentA1;
end

A=randi(60); ParentB1=Population(A,:); P1=tree_calculation(1,ParentB1,x);

B=randi(60); ParentB2=Population(B,:); P2=tree_calculation(1,ParentB2,x);

Fit_ParentB1 = sum(abs(P1-y))/1000;
 
Fit_ParentB2 = sum(abs(P2-y))/1000;

if Fit_ParentB1 > Fit_ParentB2
%ParentB=Population(index2,:); ParentB=ParentB2;
else
ParentB=ParentB1;
end

%% crossover
%create two new child
[NextGA,NextGB] = Crossover(ParentA,ParentB);
%% Mutation
%mutate two childs
NextGA = Mutation(NextGA); NextGB = Mutation(NextGB);
% compared child and parents replace with good one.
 	%calculate fitness for two child and two parents

ParentAY = tree_calculation(1,ParentA,x); ParentBY = tree_calculation(1,ParentB,x); Fit_ParentA = sum(abs(ParentAY-y))/1000; Fit_ParentB = sum(abs(ParentBY-y))/1000;

NextGAY = tree_calculation(1,NextGA,x); NextGBY = tree_calculation(1,NextGB,x); Fit_NextGA = sum(abs(NextGAY-y))/1000; Fit_NextGB = sum(abs(NextGBY-y))/1000;

% compare all of four individuals, and find best one,
 	% if child better than parents, use child, if parent better than
 	% children keep parent for next evolution.

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
%% record best fitness (MAE) for each evolution A=sort(Fit);
fitness(z) = A(1); if fitness(z) < Min
Min = fitness(z);
 



end
 
end record5(z)=Min;
DSP=[num2str(z),' : ',num2str(Min,5)]; disp(DSP);
 

function Tree=make_Tree() Tree=zeros(1,63);
% heap tree structure:
% level0 1
% level1 2,3
% level2 4-7
 
% level3 8-15
% level4 16-31
% level5 32-63
func=[11 22 33 44 55 66];
% 11 means "+", 22 means"-", 33 means", 44 means"/", 55 means “sin" 66
% means "cos", 77 means" value of x" for i = 1:63
if i== 1
R=randi(6);
Tree(i)=func(R);
elseif i== 2
R2=randi(6);
Tree(i)=func(R2);
elseif i== 3
R3=randi(6);
Tree(i)=func(R3);
elseif rand(1)<0.5 && i>3 && i<8 R4=randi(6); Tree(i)=func(R4);
elseif rand(1)<0.5 && i>7 && i<16 R5=randi(6); Tree(i)=func(R5);
elseif rand(1)<0.5 && i>15 && i<32 R4=randi(6); Tree(i)=func(R4);
else
if rand(1) < 0.5 Tree(i)= 77;
else
Num=linspace(-10,10,100000); index5=randi(100000); num=Num(index5); Tree(i)=num;
 

end end
 

end
 
end
 

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

%R= randi([2 31]); (not used)
%R= randi([2 15]);(not used)
%R= randi([2 8]);(not used) R= randi([4 15]);

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

function Tree = Mutation(Tree)
% the goal of this mutation function is that to random select a point
% between level 2 to level 4, from the sub node of start point of mutation to end,
% each node or leaf will have 50% of possibility that will replace with a
% new function or constant or X.

func = [11 22 33 44 55 66 ];
% 11 means "+", 22 means"-", 33 means", 44 means"/", 55 means “sin" 66
% means "cos", 77 means" value of x" R = randi([1 31]);
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
if rand(1) < rate Tree(R5) = 77;
else
r_range=10-15*rand(1); Tree(R5) = r_range;
end R5=R5+1;
end
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
if rand(1) < rate Tree(R4) = 77;
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
if rand(1) < rate
 
Tree(R3) = 77;
else
r_range=10-25*rand(1); Tree(R3) = r_range;
end R3=R3+1;
end

elseif R > 7 && R < 16 R1=2*R;
R2=2*R1;
for i=1:2
if rand(1) < rate num2=randi(6); Tree(R1) = func(num2);
 

end
 
end R1=R1+1;
 
for i1 = 1:4
if rand(1) < rate Tree(R2) = 77;
else
r_range=10-15*rand(1); Tree(R2) = r_range;
 

end else
 
end R2=R2+1;
 
R1=2*R;
for i = 1:2
if rand(1) < rate Tree(R1) = 77;
else
r_range=10-25*rand(1); Tree(R1) = r_range;
 


end end
 

end
 
end R1=R1+1;
 

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

