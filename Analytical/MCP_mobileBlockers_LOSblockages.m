aID = getenv('SLURM_ARRAY_TASK_ID')
if(isempty(aID))
  warning('aID is empty. Replacing it with 1.')  
  aID = '1'; %Runs only for first value of AP density when aID=1
end
rng('shuffle');

digits(10); %set precision to 10 significant digits
syms one(x); %
one(x) = kroneckerDelta(x);
syms y z w v;
syms minim(v,w) maxim(y) absol(z);
minim(v,w) = feval(symengine,'min',v,w);
maxim(y) = feval(symengine,'max',y,0);
absol(z) = abs(z);
%lij = one(max(0,abs(a*sin(theta - phi))-R))...
    %*max(0,min(r - (a*cos(theta - phi) - sqrt(R^2 - (a*sin(theta - phi))^2))...
    %,2*sqrt(R^2 - (a*sin(theta - phi))^2)));

R = 100;
V = 1; %velocity in m/s
hb = 1.8; %height blocker
hr = 1.4; %height receiver (UE)
ht = 5; %height transmitter (BS)
frac = (hb-hr)/(ht-hr);
mu = 2;

densityBS = [100,150,200,300,400,500]*10^(-6);
%densityBS = [500]*10^(-6);
densityBL_PPP = [0.01,0.1];
%densityBL_PPP = [0.1];
lambdaC = [0.01/64, 0.1/400]; %density of parent Poisson point process
% lambdaC = [0.1/400];
lambdaD = densityBL_PPP*16/pi; %density of daughter Poisson point process
rD = 1/4.0./sqrt(lambdaC); %radius of cluster disk (for daughter points)
kB = 2/pi*V.*lambdaD;

pB = 0;
N_end = 15; M_end = 15;
probBL = zeros(length(densityBS),length(densityBL_PPP),N_end,M_end);
avg_blockage_probability = zeros(length(densityBS),length(densityBL_PPP));

N_iter = 1;

for indexBS=1:length(densityBS)
     for indexBL=1:length(densityBL_PPP)  
         syms pm(m) pn(n);
         pm(m) = vpa(exp(-densityBS(indexBS)*pi*R^2)* (densityBS(indexBS)*pi*R^2)^m/factorial(m),2);
         pn(n) = vpa(exp(-lambdaC(indexBL)*pi*R^2)* (lambdaC(indexBL)*pi*R^2)^n/factorial(n),2);
         for ind_iter=1:1:N_iter
             for n=1:1:N_end
                 %     a = R*sqrt(rand(n,1));
                 %     phi = (rand(n,1)-0.5)*2*pi;
                 %     syms a [1 n];
                 %     syms phi [1 n];
                 for m=1:1:M_end
                     a = R*sqrt(rand(n,1));
                     phi = (rand(n,1)-0.5)*2*pi;
                     sprintf("n=%d m=%d",n,m)
                     syms r positive;
                     syms theta;
                     syms lij(r,theta) sumLi(r,theta) sumLi1(r,theta) sumLi2(r,theta) pBi(r,theta);
                     syms integrand(r,theta);
                     syms Rprime1(r,theta) Rprime2(r,theta) Rprime2(r,theta);
                     Rprime1(r,theta) = vpa(compose(absol,a(1)*sin(theta - phi(1)))-rD(indexBL));
                     Rprime2(r,theta) = vpa(r*frac - (a(1)*cos(theta - phi(1)) - ...
                         sqrt(compose(absol,rD(indexBL)^2 - (a(1)*sin(theta - phi(1))))^2)));
                     Rprime3(r,theta) = vpa(...
                         2*sqrt(compose(absol,rD(indexBL)^2 - (a(1)*sin(theta - phi(1)))^2)));
                     sumLi1(r,theta) = compose(compose(one,maxim),Rprime1);
                     sumLi2(r,theta) = compose(maxim,minim(Rprime2,Rprime3));
                     sumLi(r,theta) = sumLi1(r,theta)*sumLi2(r,theta);
                     if n>=2
                         for j=2:1:n
                             Rprime1(r,theta) = vpa(compose(absol,a(j)*sin(theta - phi(j)))-rD(indexBL));
                             Rprime2(r,theta) = vpa(r*frac - (a(j)*cos(theta - phi(j)) - ...
                                 sqrt(compose(absol,rD(indexBL)^2 - (a(j)*sin(theta - phi(j))))^2)));
                             Rprime3(r,theta) = vpa(...
                                 2*sqrt(compose(absol,rD(indexBL)^2 - (a(j)*sin(theta - phi(j)))^2)));
                             sumLi1(r,theta) = compose(compose(one,maxim),Rprime1);
                             sumLi2(r,theta) = compose(maxim,minim(Rprime2,Rprime3));
                             lij(r,theta) = sumLi1(r,theta)*sumLi2(r,theta);
                             sumLi = sumLi + lij;
                         end
                     end
                     pBi(r,theta) = vpa((kB(indexBL)*sumLi(r,theta))/(kB(indexBL)*sumLi(r,theta)+mu));
                     pB = (vpaintegral((vpaintegral(pBi*2*r/R^2,r,0,R))^m/2/pi,theta,-pi,pi))^m;
                     probBL(indexBS,indexBL,n,m) = probBL(indexBS,indexBL,n,m) + pB*pm(m)*pn(n);
                     %         for j=1:1:n
                     %             lij = one(max(0,abs(a(j)*sin(theta - phi(j)))-rD(indexBL)))...
                     %                 *max(0,min(r - (a(j)*cos(theta - phi(j)) - sqrt(rD(indexBL)^2 - (a(j)*sin(theta - phi(j)))^2))...
                     %                 ,2*sqrt(rD(indexBL)^2 - (a(j)*sin(theta - phi(j)))^2)));
                     
                     %         end
                 end
             end
             avg_blockage_probability(indexBS,indexBL) = sum(probBL(indexBS,indexBL,:,:),'all');
         end
         avg_blockage_probability(indexBS,indexBL) = avg_blockage_probability(indexBS,indexBL)/N_iter;
     end
end
save(strcat('AnalysisOutput','_',num2str(aID),'.mat'),'avg_blockage_probability')
%probBL = probBL/N_iter;
%prob = sum(probBL,'all')/N_iter;

