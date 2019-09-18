function [A,mu,sigma,p,smoothed] = fit_hmm(y)
T=length(y);
% Simple initial guesses for parameters - can be changed
mu=[mean(y),mean(y)]+randn(1,2)*std(y);
sigma=[std(y),std(y)];
A=[.8,.2;.2,.8];
p=.5;
iteration=2;
likelihood(1)=-999;
change_likelihood(1)=Inf;
tolerance=0.000001;
while change_likelihood(iteration-1) > tolerance
    for t=1:T % 0. probability of observing data, based on gaussian PDF
        B(t,1)=exp(-.5*((y(t)-mu(1))/sigma(1)).^2)/(sqrt(2*pi)*sigma(1));
        B(t,2)=exp(-.5*((y(t)-mu(2))/sigma(2)).^2)/(sqrt(2*pi)*sigma(2));
    end
    forward(1,:)=p.*B(1,:);
    scale(1,:)=sum(forward(1,:));
    forward(1,:)=forward(1,:)/sum(forward(1,:));
    for t=2:T % 1. probability of regimes given past data
        forward(t,:)=(forward(t-1,:)*A).*B(t,:);
        scale(t,:)=sum(forward(t,:));
        forward(t,:)=forward(t,:)/sum(forward(t,:));
    end
    backward(T,:)=B(T,:);
    backward(T,:)=backward(T,:)/sum(backward(T,:));
    for t=T-1:-1:1 % 2. probability of regime given future data
        backward(t,:)=(A*backward(t+1,:)')'.*B(t+1,:);
        backward(t,:)=backward(t,:)/sum(backward(t,:));
    end
    for t=1:T % 3-4. probability of regimes given all data
        smoothed(t,:)=forward(t,:).*backward(t,:);
        smoothed(t,:)=smoothed(t,:)/sum(smoothed(t,:));
    end
    for t=1:T-1 % 5. probability of each transition having occurred
        xi(:,:,t)=(A.*(forward(t,:)'*(backward(t+1,:).*B(t+1,:))));
        xi(:,:,t)=xi(:,:,t)/sum(sum(xi(:,:,t)));
    end
    p=smoothed(1,:);
    exp_num_transitions=sum(xi,3);
    A(1,:)=exp_num_transitions(1,:)/sum(sum(xi(1,:,:),2),3);
    A(2,:)=exp_num_transitions(2,:)/sum(sum(xi(2,:,:),2),3);
    mu(1)=(smoothed(:,1)'*y)'/sum(smoothed(:,1));
    mu(2)=(smoothed(:,2)'*y)'/sum(smoothed(:,2));
    sigma(1)=sqrt(sum(smoo  thed(:,1).*(y-mu(1)).^2)/sum(smoothed(:,1)));
    sigma(2)=sqrt(sum(smoothed(:,2).*(y-mu(2)).^2)/sum(smoothed(:,2)));
    likelihood(iteration+1)=sum(sum(log(scale)));
    change_likelihood(iteration)=abs(likelihood(iteration+1)-likelihood(iteration));
    iteration=iteration+1;
end
end
