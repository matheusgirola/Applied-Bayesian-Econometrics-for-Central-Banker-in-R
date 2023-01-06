%entering a matrix manually
X=[1 2 3 4;5 6 7 8];
Z=[2 3 4 5;7 8 9 10];
% set X(2,1)=30
X(2,1)=30;
% vertical concatenation
M=[X;Z];
%horizontal concatenation
N=[X Z];
%set the second row of N to -10
N(2,1:end)=-10;