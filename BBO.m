tic;
global succ_rate;
global feval;
global error;

for F_index=1
    
    CostFunction=@(x) test_functions(x,F_index,dim);
    [obj_val, down, up, dim, maxFE,acc_err]=test_functions_range(F_index);
    dim=2;
    VarSize=[1 dim];   % Decision Variables Matrix Size
    %      mean1=0;
    var=0;
    maxrun=1;
    MaxIt=1;  % Maximum Number of Iterations
    
    nPop=3;            % Number of Habitats (Population Size)
    succ_rate=0;
    meanerror=zeros(maxrun,1);
    meanfvalue=zeros(maxrun,1);
    meanfeval=zeros(maxrun,1);
    
    
    
    
    
    KeepRate=0.2;                   % Keep Rate
    nKeep=round(KeepRate*nPop);     % Number of Kept Habitats
    
    nNew=nPop-nKeep;                % Number of New Habitats
    
    % Migration Rates
    mu=linspace(1,0,nPop);          % Emmigration Rates
    lambda=1-mu;                    % Immigration Rates
       
    pMutation=0.7;
    
    
    
    
    %% Initialization
    for run=1:maxrun
        feval=0;
        error=0;
        
        % Empty Habitat
        habitat.Position=[];
        habitat.Cost=[];
        
        % Create Habitats Array
        pop=repmat(habitat,nPop,1);
        
        % Initialize Habitats
        for i=1:nPop
            if size(up,2)==1
                pop(i).Position=unifrnd(down,up,VarSize);
            end
            if size(up,2)>1
                for j=1:dim
                    high=up(j);low=down(j);
                    pop(i).Position(j)=unifrnd(high,low,[1 1]);
                end
            end
            
            pop(i).Cost=test_functions(pop(i).Position,F_index,dim);
            feval=feval+1;
           
        end
        
        for i=1:nPop
            
                disp(['Island',num2str(i),': ',num2str(pop(i).Position),'    cost: ',num2str(pop(i).Cost)]);
        end
        % Sort Population
        [~, SortOrder]=sort([pop.Cost]);
        pop=pop(SortOrder);
        disp(['Sorted ']);
        for i=1:nPop
            
                disp(['Island',num2str(i),': ',num2str(pop(i).Position),'    cost: ',num2str(pop(i).Cost)]);
        end
        
        % Best Solution Ever Found
        BestSol=pop(1);
        
        % Array to Hold Best Costs
        BestCost=zeros(MaxIt,1);
        
        
        
        % BBO Main Loop
        % while FE<maxFE
        
        for it=1:MaxIt
            
            newpop=pop;
            for i=1:nPop
                for k=1:dim
                    prob=rand;
                    % Migration
                    if prob<=lambda(i)
                        % Emmigration Probabilities
                        disp(['Replace ',num2str(i),',',num2str(k)]);
                        EP=mu;
                        EP(i)=0;
                        EP=EP/sum(EP);
                        
                        % Select Source Habitat
                        j=RouletteWheelSelection(EP);
                        
                        disp(['From ',num2str(j),',',num2str(k)]);
%                         disp(['pop(j).Position(j) ',num2str(pop(j).Position(k))]);
%                         disp(['pop(i).Position(j) ',num2str(pop(i).Position(k))]);
%                         disp(['pop(j).Position(k)-pop(i).Position(k) ',num2str(pop(j).Position(k)-pop(i).Position(k))]);
                        
                        newpop(i).Position(k)=pop(j).Position(k); 
                            
%                         disp(['trumig ',num2str(i),',',num2str(k),'tru:',num2str(newpop(i).Position(k))]);
                    end
                    
                    %% Mutation
%                     if size(up,2)==1
%                         sigma=0.02*(up-down);
%                         
%                     end
%                     if size(up,2)>1
%                         sigma=0.02*(up(k)-down(k));
%                     end
                    temprand2=rand;
                    pMutation;
                    if temprand2<=pMutation
                        if size(up,2)==1
                            newpop(i).Position=unifrnd(down,up,VarSize);
                        end
                        if size(up,2)>1
                            for j=1:dim
                                high=up(k);low=down(k);
                                newpop(i).Position(k)=unifrnd(high,low,[1 1]);
                            end
                        end
                       disp(['pmut true: ',num2str(i),',',num2str(k),'pmut:',num2str(newpop(i).Position(k))]);
                    end
                end
            
                % Apply Lower and Upper Bound Limits
                newpop(i).Position = max(newpop(i).Position, down);
                newpop(i).Position = min(newpop(i).Position, up);
                
                % Evaluation
                newpop(i).Cost=test_functions(newpop(i).Position,F_index,dim);
                feval=feval+1;
                
            end
            disp(['new pop ']);
            for i=1:nPop
            
                disp(['Island',num2str(i),': ',num2str(pop(i).Position),'    cost: ',num2str(pop(i).Cost)]);
            end
            % Sort New Population
            [~, SortOrder]=sort([newpop.Cost]);
            newpop=newpop(SortOrder);
            disp(['Sorted new pop ']);
            for i=1:nPop
            
                disp(['Island',num2str(i),': ',num2str(pop(i).Position),'    cost: ',num2str(pop(i).Cost)]);
            end
            
            % Select Next Iteration Population
            pop=[pop(1:nKeep)
                newpop(1:nNew)];
            disp(['Combined pop for next gen ']);
            for i=1:nPop
            
                disp(['Island',num2str(i),': ',num2str(pop(i).Position),'    cost: ',num2str(pop(i).Cost)]);
            end
            % Sort Population
            [~, SortOrder]=sort([pop.Cost]);
            pop=pop(SortOrder);
            disp(['Sorted pop for next gen ']);
            for i=1:nPop
            
                disp(['Island',num2str(i),': ',num2str(pop(i).Position),'    cost: ',num2str(pop(i).Cost)]);
            end
            % Update Best Solution Ever Found
            BestSol=pop(1);
                      
            %     Store Best Cost Ever Found
            BestCost(it)=BestSol.Cost;
          
            error=abs(BestSol.Cost-obj_val); %error
            if ((abs( BestSol.Cost-obj_val)<=acc_err)||( BestSol.Cost<obj_val)) % succ_rate
                succ_rate= succ_rate+1;
                disp([' : succ_rate ' num2str(succ_rate) ' : iter' num2str(it)  ' : run ' num2str(run) ]);
                break;
            end
                      
            
            if feval>= maxFE
                
                break;
            end
            
            
        end
        
        meanerror(run)=error;
        meanfeval(run)=feval;
        meanfvalue(run)=BestSol.Cost;
        
   end    % end of run loop
    
    ferror=mean(meanerror);
    fvalue=mean(meanfvalue);
    ffeval=mean(meanfeval);
    sd=std(meanfvalue);
    
    
    
    
    
    % figure;
    % %plot(BestCost,'LineWidth',2);
    % semilogy(BestCost,'LineWidth',2);
    % xlabel('Iteration');
    % ylabel('Best Cost');
    % grid on;
    disp(['F_index=' num2str(F_index)  'Total Run=' num2str(maxrun), ':succ_rate=' num2str(succ_rate) 'Mean Feval=' num2str(ffeval)   ':Mean Fun Val=' num2str(fvalue) ':Mean Error=' num2str(ferror) ':Std=' num2str(sd) ]);
    
    
    filename='bbo.xlsx';
    F_index = [F_index]; dim=[dim];Total_run=[run]; succ_rate=[succ_rate]; Mean_Feval=[ffeval]; Mean_Fun_Val=[fvalue]; Mean_Error=[ferror];Std=[sd];acc_err=[acc_err];
    fileExist = exist(filename,'file');
    if fileExist==0
        header = {'F_index','dim', 'run_num ','succ_rate' , 'Mean Feval','Mean Fun Val','Mean Error','Std','acc_err'};
        xlswrite(filename,header);
    end
    [~,~,input] = xlsread(filename); % Read in your xls file to a cell array (input)
    new_data = {F_index,dim,Total_run ,succ_rate ,Mean_Feval,Mean_Fun_Val,Mean_Error,Std,acc_err}; % This is a cell array of the new line you want to add
    output = cat(1,input,new_data); % Concatinate your new data to the bottom of input
    xlswrite(filename,output); % Write to the new excel file.
    
end
