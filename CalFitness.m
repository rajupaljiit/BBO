
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 
function fFitness=CalFitness(fitness) 
fFitness=zeros(size(fitness));
ind=find(fitness>=0);
fFitness(ind)=1./(fitness(ind)+1);
ind=find(fitness<0);
fFitness(ind)=1+abs(fitness(ind));

%function fFitness=CalFitness(fitness)
%fFitness=zeros(size(fObjV));