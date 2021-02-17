clearvars close all clc
%% initialization Dataset=load('SR_div_1000.txt'); x=Dataset(:,1)';
y=Dataset(:,2)';
MAE = 100; %min_error %min mean error fitness=zeros(1,100000);

%% start tree. random generate a heap tree Tree=make_Tree();
HCTree=Tree;
%% bring initial tree into modification. for i=1:100000
 	% in HC_Tree function, random a start points between root (1) to last
 	% point (63), from start point to its sub node, all of function nodes leaf’s have 30%
 	% possibility switch, function will be switched to another function,
 	% constant and x will switch a new random constant or x.
Tree=HC_Tree(HCTree);
% calculate the Y Y0=tree_calculation(1,Tree,x);
% calculate mean absolute value (MAE) fitness(i)= sum(abs(Y0-y))/1000;
 	% if fitness is better than the previous one, this result will be new MAE if fitness(i)<MAE
HCTree=Tree; % update tree. MAE=fitness(i);
 


end
 
end Fitness(i)=MAE; disp(MAE);
 

%generate a tree. function Tree=make_Tree() Tree=zeros(1,63);
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
if i== 1 %root level R=randi(6); Tree(i)=func(R);
elseif i== 2 %1 level R2=randi(6); Tree(i)=func(R2);
elseif i== 3 %1 level R3=randi(6); Tree(i)=func(R3);
elseif rand(1)<0.5 && i>3 && i<8 %2 level R4=randi(6);
Tree(i)=func(R4);
elseif rand(1)<0.5 && i>7 && i<16 %3 level R5=randi(6);
Tree(i)=func(R5);
elseif rand(1)<0.5 && i>15 && i<32 %4 level R4=randi(6);
Tree(i)=func(R4);
else %if rand there are no node, set up as leaf and leaf will be whether constant or x if rand(1) < 0.5
Tree(i)= 77; else
Num=linspace(-10,10,100000); index5=randi(100000); num=Num(index5); Tree(i)=num;
 

end end
 

end
 
end
 

% Recursion method to calculate the Y from 1 to 1000 x
% Recursion will search the tree from top to bottom to identify what each node and leaf,
% then will start to calculate from bottom to top, the final results will be Y(i) function Y0=tree_calculation(i,Tree,x)
block=Tree(i);
if block == 11 %%if i=11, it bring 2*i and 2i+1 node to add each other, Y0=tree_calculation(2*i,Tree,x) + tree_calculation(2*i+1,Tree,x);
elseif block == 22
Y0=tree_calculation(2*i,Tree,x) - tree_calculation(2*i+1,Tree,x); elseif block == 33
Y0=tree_calculation(2*i,Tree,x) .* tree_calculation(2*i+1,Tree,x); elseif block == 44
Y0=tree_calculation(2*i,Tree,x) ./ tree_calculation(2*i+1,Tree,x); elseif block == 55
Y0=sin(tree_calculation(2*i,Tree,x)); elseif block == 66
Y0=cos(tree_calculation(2*i,Tree,x)); elseif block == 77
Y0= x;
else
Y0=Tree(i);
end end

function Tree = HC_Tree(Tree)
% the goal of this mutation function is that to random select a point
% between level 2 to level 4, from the sub node of start point of mutation to end,
% each node or leaf will have 50% of possibility that will replace with a
% new function or constant or X.

func = [11 22 33 44 55 66 ];
% 11 means "+", 22 means"-", 33 means", 44 means"/", 55 means “sin" 66
% means "cos", 77 means" value of x"
 
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
end R1=R1+1;
end
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
end R1=R1+1;
end
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
end R4=R4+1;
end

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
if rand(1) < rate Tree(R3) = 77;
else
r_range=10-15*rand(1); Tree(R3) = r_range;
 

end
 
end R3=R3+1;
 

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
r_range=10-15*rand(1);
 
Tree(R1) = r_range;
end
 

end end
 

end
 
R1=R1+1;
