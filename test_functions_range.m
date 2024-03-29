% GSA code v1.1.
% Generated by Esmat Rashedi, 2010. 
% "	E. Rashedi, H. Nezamabadi-pour adim S. Saryazdi,
%�GSA: A Gravitational Search Algorithm�, Information sciences, vol. 179,
%no. 13, pp. 2232-2248, 2009."
%
% This function gives boudimaries adim dimension of search space for test functions.

function [obj_val, down, up, dim, maxFE,acc_err]=test_functions_range(F_index)
maxFE=200000; % maximum number of function evaluatiodim
% Tol = 1e-12;
% Benchmark function
if F_index==1      % Parabola (Sphere)
    down=0;
    up=pi/2;
    dim=3; 
   obj_val=0;
   acc_err=1.0e-1;
   
elseif F_index==2       %De Jong's f4
     down=-5.12;
    up=5.12;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==3        %Griewank
   down=-600;
    up=600;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==4        %Rosenbrock
  down=-100;
    up=100;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-2;
    
elseif F_index==5       %Rastrigin. Minimum value 0. Solution (0,0 ...0)
   down=-5.12;
    up=5.12;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==6       % Ackley
    down=-30;
    up=30;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==7       %Alpine
  down=-10;
    up=10;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==8       %Michalewicz
  down=0;
    up=pi;
    dim=10; 
   obj_val = -9.66015;
   acc_err=1.0e-5;
    
elseif F_index==9       %Cosine Mixture [-1,1] f(0,0,...0)=-D*0.1
   down=-1;
    up=1;
    dim=30; 
   obj_val = -dim*0.1;
   acc_err=1.0e-5;
    
elseif F_index==10      %Exponential [-1,1] f(0,0,...0)= -1
    down=-1;
    up=1;
    dim=30; 
   obj_val = -1;
   acc_err=1.0e-5;
   
elseif F_index==11      %Zakharov's [-5.12,5.12] f(0,0,0,...,0)=0
 down=-5.12;
    up=5.12;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-2;
  
elseif F_index==12      %Cigar [-10,10]  f(0,0,0,...,0)=0
   down=-10;
    up=10;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==13      %brown3 [-1,4]  f(0,0,0.....,0)=0
   down=-1;
    up=4;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
   
elseif F_index==14      %/Schewel prob 3
   down=-10;
    up=10;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==15      %Salomon Problem (SAL)
  down=-100;
    up=100;
    dim=30; 
   obj_val = 0;
   acc_err=2.0e-1;
    
elseif F_index==16      %Axis parallel hyperellipsoid
   down=-5.12;
    up=5.12;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==17      %Pathological
    down=-100;
    up=100;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-1;
    
elseif F_index==18      %Sum of different powers
   down=-1;
    up=1;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==19      %step function [-100, 100] f(-0.5<=x<=0.5)=0
    down=-100;
    up=100;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==20      %Quartic function, i.e., noise [-1.28, 1.28] f(0000..00)=0
  down=-1.28;
    up=1.28;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==21      %Inverted cosine wave function (Masters) //Inverted cosine wave function (Masters) [-5, 5] f(000..0)=-D+1
    down=-5;
    up=5;
    dim=10; 
   obj_val = -dim+1;
   acc_err=1.0e-5;
   
elseif F_index==22      %Neumaier 3 Problem (NF3) (Neumaier, 2003b)
    down=-900;
    up=900;
    dim=10; 
   obj_val = -(dim*(dim+4)*(dim-1))/6.0;
   acc_err=1.0e-5;
    
elseif F_index==23      %Rotated hyper-ellipsoid function
  down=-65.536;
    up=65.536;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==24      %Levi montalvo 1
     down=-10;
    up=10;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==25        %Levi montalvo 2 
   down=-5;
    up=5;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
   
elseif F_index==26      %Ellipsoidal Ellipsoidal [-D,D] f(1,2,3,...,D)=0
    down=-30;
    up=30;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-5;
    
elseif F_index==27      %Beale function [-4.5,4.5] f(3, 0.5)=0
     down=-4.5;
    up=4.5;
    dim=2; 
   obj_val = 0;
   acc_err=1.0e-5;
   
elseif F_index==28      %Colville function [-10,10] f(1111)=0
     down=-10;
    up=10;
    dim=4; 
   obj_val = 0;
   acc_err=1.0e-3;
   
elseif F_index==29      %Branins�s function [-5,10][0,15] f(-pi, 12.275)=0.3979
    down=[-5, 0];
    up=[10, 15];
    dim=2; 
   obj_val = 0.3979;
   acc_err=1.0e-4;
    
elseif F_index==30      %Kowalik function [-5,5] f(0.192833, 0.190836, 0.123117, 0.135766)=0.000307486
   down=-5;
    up=5;
    dim=4; 
   obj_val = 0.000307486;
   acc_err=1.0e-4;  
    
%%%%%%// Minibench Mark Problems: 4-Problem Taken from clerc website%%%%%%
    
elseif F_index==31      %2D Tripod function [-100,100] f(0, 50)=-50
    down=-100;
    up=100;
    dim=2; 
   obj_val = 0;
   acc_err=1.0e-6;
   
elseif F_index==32     %Shifted CEC 2005 Rosenbrock F6  [-100, 100], solution point is O + (1; : : : ; 1) where f = 390.
   down=-100;
    up=100;
    dim=10; 
   obj_val = 390;
   acc_err=1.0e-1;
    
elseif F_index==33      %Shifted Parabola/Sphere (CEC 2005 benchmark)	x?[-100,100] , Global optimum: x* = offset  f(x) = f_bias = - 450
    down=-100;
    up=100;
    dim=10; 
   obj_val = -450;
   acc_err=1.0e-5; 
	    
elseif F_index==34      %Shifted CEC 2005  Rastrigin  x?[-5,5] , Global optimum x* = offset , f(x*) = f_bias = - 330
    down=-5;
    up=5;
    dim=10; 
   obj_val = -330;
   acc_err=1.0e-2;
    
elseif F_index==35      %Shifted CEC 2005 Schwefel [-100,100], Global optimum x* = offset , f(x*) = f_bias = - 450
    down=-100;
    up=100;
    dim=10; 
   obj_val = -450;
   acc_err=1.0e-5;
          
elseif F_index==36        %Shifted CEC 2005 Griewank. WARNING: in the CEC 2005 benchmark it is rotated
     down=-600;
    up=600;
    dim=10; 
   obj_val = -180;
   acc_err=1.0e-5;
    
elseif F_index==37      %Shifted Ackley (CEC 2005) [-32,32], Global optimum x* = offset , f(x*) = f_bias = - 140
     down=-32;
    up=32;
    dim=10; 
   obj_val = -140;
   acc_err=1.0e-5;
   
% elseif F_index==38        %Compression spring [1,...,70],[0.6,3],[0.207,0.5] f(7; 1:386599591; 0:292) = 2:6254214578.
%     
%     down=[1, 0.6, 0.207];
%     up=[70, 3, 0.5];
%     dim=3; 
%    obj_val = 2.6254214578;
%    acc_err=1.0e-10;

elseif F_index==39       % Gear Train Problem
    down=12;
    up=60;
    dim=4; 
   obj_val = 2.7e-12;
   acc_err=1.0e-15;

elseif F_index==40      %Goldstein-Price function  [-2,2] f(0, -1) = 3.
     down=-2;
    up=2;
    dim=2; 
   obj_val = 3;
   acc_err=1.0e-14;
   
elseif F_index==41    %Six-hump camel back function  [-5 < x < 5]
     down=[-5, -5];
    up=[5, 5];
    dim=2; 
   obj_val = -1.0316;
   acc_err=1.0e-5;
   
elseif F_index==42      %Easom's function  -10<=x(i)<=10, i=1:2. f(x1,x2)=-1; (x1,x2)=(pi,pi).
     down=-100;
    up=100;
    dim=2; 
   obj_val = -1;
   acc_err=1.0e-13;
    
elseif F_index==43      %Dekkers and Aarts Problem (DA)  -20<=x(i)<=20, i=1:2. f(0,15),f(0,-15)=-24777; 
    down=-20;
    up=20;
    dim=2; 
   obj_val = -2477;
   acc_err=0.5;

elseif F_index==44      %Hosaki Problem (HSK)   0<=x1<=5, 0<=x2<=6,f(4,2)=-2.3458;
   down=[0, 0];
    up=[5, 6];
    dim=2; 
   obj_val = -2.3458;
   acc_err=1.0e-5;

elseif F_index==45      %McCormick Problem (MC)   -1.5<=x1<=4, -3<=x2<=3,f(4,2)=-2.3458; 
    down=[-1.5, -3];
    up=[4, 3];
    dim=2; 
   obj_val = -1.9133;
   acc_err=1.0e-4;
  
elseif F_index==46      %Meyer and Roth Problem (MR) (Wolfe, 1978)  -10<=x<=10, f(�3.13; 15.16; 0.78)=0.4*10^(-4);    
     down=-10;
    up=10;
    dim=3; 
   obj_val = 0.4e-4;
   acc_err=1.95e-3;

elseif F_index==47      %Shubert Problem (SBT)  -10<=x<=10, f((7.0835, 4.8580)=-186.7309;
    down=-10;
    up=10;
    dim=2; 
   obj_val = -186.7309;
   acc_err=1.0e-5;
    
%elseif F_index==48      %Sinusoidal Problem (SIN) 0<=x<=180, A= 2.5; B =5; z = 30. f(90 + z; 90 + z; . . . ; 90 + z)= -(A+1)
      down=-0;
    up=180;
    dim=10; 
   obj_val = -3.5;
   acc_err=1.0e-2;
   
elseif F_index==49      % Moved axis parallel hyper-ellipsoid function [-5.12, 5.12] f(x)=0; x(i)= 5*i, i=1:D.
    down=-5.12;
    up=5.12;
    dim=30; 
   obj_val = 0;
   acc_err=1.0e-15;
    
elseif F_index==50      % Pressure vessel (confinement method) 	// (1.125, 0.625, 55.8592, 57.7315) => 7197.729
			% If no granularity => min = 6059.7143
			% We are using function without granularity 
     down=[1.125, 0.625, 1.0e-8, 1.0e-8];
    up=[12.5, 12.5, 240, 240];
    dim=4; 
   obj_val = 7197.729;
   acc_err=1.0e-5;
    
% elseif F_index==51      %lennard_jones (-2,2) for 5 atoms => -9.103852
%   down=-2;
%     up=2;
%     dim=15; 
%    obj_val = -9.103852;
%    acc_err=1.0e-4;

% elseif F_index==52      %Parameter Estimation for Frequency-Modulated (FM) Sound Waves f(x)=0;
%   down=-6.4;
%     up=6.35;
%     dim=6; 
%    obj_val = 0;
%    acc_err=1.0e-2;
   
elseif F_index==53      %temp
    down=-5;
    up=5;
    dim=2; 
   obj_val = 0;
   acc_err=1.0e-5;
end
