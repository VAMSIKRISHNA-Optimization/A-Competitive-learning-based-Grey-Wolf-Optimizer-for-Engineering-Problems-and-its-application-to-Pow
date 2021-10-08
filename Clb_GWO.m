


function [Alpha_score,Alpha_pos,Convergence_curve]=Clb_GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems
fitness_new=inf;
%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);


for i=1:size(Positions,1)
    Fit(i) = fobj(Positions(i,:));
end

pBestScore = Fit;
pBest = Positions;




Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter
nn=0.5;
Fail=0;
pAss=0;

% Main loop
while l<Max_iter
    
    Pop1=round(0.9*SearchAgents_no);
    Pop2=round(0.1*SearchAgents_no);
    
    for i=1:size(Positions,1)
        
       
        
        % Calculate objective function for each search agent
        fitness = Fit(i);
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<Beta_score 
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
          [~,ind] = max(fitness);
               Worst = Positions(ind,:);
    end
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
       f1=rand;
       f2=rand;
       f3=rand;
            r1=rand();
            r2=rand();
        for j=1:size(Positions,2)     
            
            if i<Pop1
            
            
         
             Sel1=rand;
             if Sel1>0.5   
            %Standard GWO scheme

            
            A1=2*a*r1-a; 
            C1=2*r2; 
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); 
            X1=Alpha_pos(j)-A1*D_alpha; 
                       
       
            
            A2=2*a*r1-a; 
            C2=2*r2; 
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); 
            X2=Beta_pos(j)-A2*D_beta; 
            A3=2*a*r1-a; 
            C3=2*r2; 
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j));  
            X3=Delta_pos(j)-A3*D_delta;      
            
            OldPositions(i,j)=(X1+X2+X3)/3;
             else
           % Non-linear exploration scheme
           R1=randi(Pop1,1,1);
           R2=randi(Pop1,1,1);
           R3=randi(Pop1,1,1);
           

S=exp(-0.005*l); %Same as Phi in the manuscript

            Y1=S.*Positions(R1,j)-f1.*((Alpha_pos(j)-Worst(j))-Positions(i,j));
            Y2=S.*Positions(R2,j)-f2.*((Beta_pos(j)-Worst(j))-Positions(i,j));
            Y3=S.*Positions(R3,j)-f3.*((Delta_pos(j)-Worst(j))-Positions(i,j));
            
            OldPositions(i,j)=(Y1+Y2+Y3)/3;
            
             end
            else
                         Sel2=rand;
                         if Sel2<0.5
                                   R1=randi(Pop1,1,1);
                                   R2=randi(Pop1,1,1);
                        OldPositions(i,j)=Positions(i,j)+(1+rand).*(Positions(R1,j)-Positions(R2,j));   %Re-intitialization of 10% of population
                         else
                        OldPositions(i,:)=lb+rand*(ub-lb); 
                         end
                    end
        
    end
   
        
     
        
        
        
           % Return back the search agents that go beyond the boundaries of the search space
           OldPositions(i,:)=min(ub, OldPositions(i,:));
           OldPositions(i,:)=max(lb, OldPositions(i,:));
        
  

               Fit_GWO(i)=fobj(Positions(i,:));

  if  Fit_GWO(i)<Alpha_score
            Alpha_score=Fit_GWO(i); % Update alpha
            Alpha_pos=OldPositions(i,:);
        end
        
        if Fit_GWO(i)>Alpha_score && Fit_GWO(i)<Beta_score 
            Beta_score=Fit_GWO(i); % Update beta
            Beta_pos=OldPositions(i,:);
        end
        
        if Fit_GWO(i)>Alpha_score && Fit_GWO(i)>Beta_score && Fit_GWO(i)<Delta_score 
            Delta_score=Fit_GWO(i); % Update delta
            Delta_pos=OldPositions(i,:);
        end


    end
    
    % Competetive learning phase
    
    
    for i=1:SearchAgents_no
        
         if i<Pop2
              p = randi([1 SearchAgents_no],1,1);      % Selection of random parter
        
        %% Ensuring that the current member is not the partner
        while i == p
            p = randi([1 SearchAgents_no],1,1);  % Selection of random parter
        end
             if(Fit_GWO(i)<Fit_GWO(p))
                        Xnew(i,:)=OldPositions(p,:)+rand.*(Positions(i,:)-Positions(p,:))+rand.*(Alpha_pos-Positions(p,:));

                        else
                            Xnew(i,:)=Positions(i,:)-rand.*(Alpha_pos-Positions(p,:));
             end
         else



        Candidates = [1:i-1 i+1:SearchAgents_no];            % Ensuring that the current member is not the partner
        idx = Candidates(randperm(SearchAgents_no-1,5));     % Selection of three random parters
        
        R1 =Positions(idx(1),:);                       % Assigning randomly selected solution 1
        R2 = Positions(idx(2),:);                      % Assigning randomly selected solution 2
        R3 = OldPositions(idx(3),:);                % Assigning randomly selected solution 3 from the hunting scheme
         R4 = Positions(idx(4),:);    
         R5 = Positions(idx(5),:);    
  
        
        
        

          selsc=rand();

          if selsc>0.75
              ne=rand;
              if ne>0.75
                     Xnew(i,:)= Alpha_pos+rand*(R1-R2)+rand*(R3-R4);  
              else

            Xnew(i,:)= Alpha_pos+a*(R1-R2);  %Linear learning
              end
           
            
            else

            
            if nn>0.5
            Xnew(i,:)=OldPositions(i,:)+rand(1,dim).*(R3-Alpha_pos)-rand(1,dim).*(OldPositions(i,:)+(Beta_pos+Delta_pos));     % Adaptive learning

            else
                
             Xnew(i,:)=R3+rand.*(R1-R2)+rand.*(R4-R5);
            end

        
          end
         end
       Xnew(i,:)=min(ub,   Xnew(i,:));
       Xnew(i,:)=max(lb,   Xnew(i,:));
        
        
        

         Fit_Clb(i)=fobj(Xnew(i,:));
    end

     
 %% Adaptive mechanism
    tmp = Fit_GWO < Fit_Clb;                         
    if tmp==1
        Fail=Fail+1;
        if Fail==20
       
        nn=rand;
        Fail=0;
        end
    end
    
    
        if tmp==0
        pAss=pAss+1;
        if pAss==10
       pAss=Fail;
        nn=abs(1-nn);
        end
    end
    tmp_rep = repmat(tmp',1,dim);
    
    tmpFit = tmp .* Fit_GWO + (1-tmp) .* Fit_Clb;
    tmpPositions = tmp_rep .* OldPositions + (1-tmp_rep) .* Xnew;
    
    %% Updating
    tmp = pBestScore <= tmpFit;                             
    tmp_rep = repmat(tmp',1,dim);
    
    pBestScore = tmp .* pBestScore + (1-tmp) .* tmpFit;
    pBest = tmp_rep .* pBest + (1-tmp_rep) .* tmpPositions;
    
    Fit = pBestScore;
    Positions = pBest;
    
   
       l=l+1; 



disp([ l abs(Alpha_score)])
Convergence_curve(l)=Alpha_score;

end
end




