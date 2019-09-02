function prob = CalProb(fFitness(i))
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

prob=(0.9.*fFitness(i)./max(fFitness))+0.1;
%disp('prob');disp(prob);
end

