function [Destination_fitness,bestPositions,Convergence_curve]=RCLSMAOA(N,Max_iter,lb,ub,dim,fobj) 
bestPositions=zeros(1,dim);         % initialize position
Destination_fitness=inf;            % For maximization issues, change it to - inf
AllFitness = inf*ones(N,1);         % Record the fitness of all slime molds
weight = ones(N,dim);               % Suitable weight for each type of slime mold
X=initialization(N,dim,ub,lb);      % Initialize Random Solution Set
Xnew=X;
Xnew_fobj=zeros(1,N);
Xc=(ub+lb)/2;
tri=zeros(1,N);
Convergence_curve=zeros(1,Max_iter);
it=1;                               % Number of iterations
lb=ones(1,dim).*lb;                 % lower boundary 
ub=ones(1,dim).*ub;                 % upper boundary
z=0.03;                             % parameter
Alpha=5;
Mu=0.499;                          % Main loop
while  it <= Max_iter
    limit=log(it);
    for i=1:N
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        AllFitness(i) = fobj(X(i,:));
    end
    
    %Select the best and worst fitness values
    [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
    worstFitness = SmellOrder(N);
    bestFitness = SmellOrder(1);
    
    % SMA
    S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero
    % AOA
    MOP=1-((it)^(1/Alpha)/(Max_iter)^(1/Alpha));   % Probability Ratio
    
    
    % Update weights
    for i=1:N
        for j=1:dim
            if i<=(N/2)
                weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            else
                weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
            end
        end
    end
    
    %update the best fitness value and best position
    if bestFitness < Destination_fitness
        bestPositions=X(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    
    a = atanh(-(it/Max_iter)+1);
    
    % Update the Position of search agents
    for i=1:N
        if rand<z     %
            Xnew(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(AllFitness(i)-Destination_fitness));
            vb = unifrnd(-a,a,1,dim);
            r = rand();
            A = randi([1,N]);  % two positions randomly selected from population
            B = randi([1,N]);
            if r<p    %Eq.(2.1)
                Xnew(i,:) = bestPositions+ vb.*(weight(i,:).*X(A,:)-X(B,:));
            else
                r1=rand();
                if r1>0.5
                    Xnew(i,:)=bestPositions./(MOP+eps).*((ub-lb).*Mu+lb);
                else
                    Xnew(i,:)=bestPositions.*MOP.*((ub-lb).*Mu+lb);
                end
            end
            Xnew_fobj(1,i)=fobj(Xnew(i,:));
        end
        if Xnew_fobj(1,i) < AllFitness(i)
            X(i,:)=Xnew(i,:);
            AllFitness(i)= Xnew_fobj(1,i);
            if AllFitness(i)< bestFitness
                bestPositions=X(i,:);
                bestFitness=AllFitness(i);
            end
        end
        
    end
    
    for i=1:N    %随机中心解
        xl = ceil(N*rand);
        while xl==i
            xl = ceil(N*rand);
        end
        xr = ceil(N*rand);
        while xr==i
            xr = ceil(N*rand);
        end
        Xc_rand=rand;
        if Xc_rand<0.5
            Xnew(i,:)=Xc+(X(xr,:)-Xc).*rand();
        else
            Xnew(i,:)=Xc+(Xc-X(xl,:)).*rand();
        end
        
    end
    
    
    
    
    %Mutation Strategy
    for i=1:N
        v1 = X(i,:);
        v2 = X(i,:);
        v3 = X(i,:);
        for j = 1:dim
            if rand()<0.1 || j ==randi(dim)
                v1(1,j) = X(randi(N),j)+1.0*(X(randi(N),j)-X(randi(N),j));
            end
            if rand()<0.2 || j ==randi(dim)
                v2(1,j) = X(randi(N),j)+0.8*((X(randi(N),j)-X(randi(N),j)))+0.8*((X(randi(N),j)-X(randi(N),j)));
            end
            if rand()<0.9 || j ==randi(dim)
                v3(1,j) = X(i,j)+rand()*(X(randi(N),j)-X(i,j))+(X(randi(N),j)-X(randi(N),j));
            end
        end
        a = fobj(v1(1,:));
        b = fobj(v2(1,:));
        c = fobj(v3(1,:));
        
        if a<b
            if a<c
                Xnew(i,:)=v1(1,:);
            else
                Xnew(i,:)=v3(1,:);
            end
            
        else
            if c<b
                Xnew(i,:)=v2(1,:);
            else
                Xnew(i,:)=v3(1,:);
            end
        end
    end
    
    %Restart strategy
    for i=1:N
        if tri(i)<limit
            if fobj(Xnew(i,:))<AllFitness(i)
                tri(i)=0;
            else
                tri(i)=tri(i)+1;
            end
        else
            t1 = zeros(1,dim); t2 = zeros(1,dim);
            t1(1,:) = (ub-lb)*rand()+lb;
            t2(1,:) = (ub+lb)*rand()-Xnew(i,:);
            Flag4ub=t2(1,:)>ub;
            Flag4lb=t2(1,:)<lb;
            t2(1,:)=(t2(1,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;%拉回至边界位置
            if fobj(t1) < fobj(t2)
                Xnew(i,:) = t1(1,:);
            else
                Xnew(i,:) = t2(1,:);
            end
            tri(i)=0;
            
        end
    end
    Convergence_curve(it)=Destination_fitness;
    it=it+1;
end
end