% This is a general program to perform a Rainfall Cascade Evaluation
% The program... 

clear
close all

% load data, assumed to be in textfile with one column and nr rows, where nr is a number 2^n

%file='Mumbai_colaba_15min_fix.dat';
%file='B.txt'; %we have used 45 min values as the data to be broken down is in 24 hours
%tres=45; %time resolution of the data is 45 min
tres=900; %time resolution of the data in seconds
load A
F11=A(:,11);
save('F11.mat','F11');
load F11.mat
% set the number of "cascade steps" to be evaluated
% starting from the highest resolution

nstep=7; %cascade steps-5 because the data needs to be broken down in 45 min instead of 15 min because
% in 5 cascade steps it will brake to 45 min and 6 step it goes to 22.5 min which leads it to crash the script
boxs=[1;2;4;8;16;32;64;128;256;512;1024];
nvc=3; %volume classes 
boxt=boxs*tres;
for cs=1:nstep
    cst(cs)=mean(boxt(cs:cs+1));
end
inname=F11;
% inname=['D:\test fractal\',file];
data(:,1)=inname;

[nr,nc] = size(data); % nc - columns (1), nr - rows

% make a matrix with successively 2-by-2 aggregated values
% it has nr rows in the first column, then successively halved nr, 
% and nstep+1 columns

for cs=1:nstep
    npar=nr/(2^cs); %number of "parents"
    for ap=1:npar %ap=aktuell "parent"
        ao=ap*2; %ap=aktuell "offspring"
        data(ap,cs+1)=data(ao-1,cs)+data(ao,cs);
    end
end

% start by creating matrices with "parent information"

for cs=1:nstep %nstep
    npar=nr/(2^cs);
    parn(1:4)=0;
    clear parinfo;
    for ap=2:npar-1 %skip first and last as the position type is unknown
        if (data(ap,cs+1)>0)            
            if ((data(ap-1,cs+1)==0) && (data(ap+1,cs+1)==0))
                parn(1)=parn(1)+1;
                parval(parn(1),1)=data(ap,cs+1); %matrix with "parent values"
                parinfo(ap,1)=1; %column 1 "parent types": 1 - isolated, 2 - starting, 3 - enclosed, 4 - stopping
            elseif ((data(ap-1,cs+1)==0) && (data(ap+1,cs+1)>0))
                parn(2)=parn(2)+1; 
                parval(parn(2),2)=data(ap,cs+1);    
                parinfo(ap,1)=2;                 
            elseif ((data(ap-1,cs+1)>0) && (data(ap+1,cs+1)>0))
                parn(3)=parn(3)+1;
                parval(parn(3),3)=data(ap,cs+1);    
                parinfo(ap,1)=3;                    
            else
                parn(4)=parn(4)+1;
                parval(parn(4),4)=data(ap,cs+1);    
                parinfo(ap,1)=4;                    
            end
            ao=ap*2;
            if (data(ao-1,cs)==0)
               parinfo(ap,3)=1; %column 3 "division types": 1 - 0/1, 2 - 1/0, 3 - x/x
            elseif (data(ao,cs)==0)
               parinfo(ap,3)=2;
            else
               parinfo( ap,3)=3;     
               parinfo(ap,4)=data(ao-1,cs)/data(ap,cs+1); %column 4 "x/x-weight"
            end
        end
    end

    parn;

% based on four volume classes, make volume limits (this is for 4 volume
% clases)
                
%   lim(1,1:3)=prctile(parval(1:parn(1),1),[25 50 75]);  
%   lim(2,1:3)=prctile(parval(1:parn(2),2),[25 50 75]);  
%   lim(3,1:3)=prctile(parval(1:parn(3),3),[25 50 75]);  
%   lim(4,1:3)=prctile(parval(1:parn(4),4),[25 50 75]);     
    
%   lim

% add volume class in parinfo, column 2: 1 - smallest,...,5 - highest
    
%   for pc=1:4
%       parpos=find(parinfo(:,1)==pc); %positions of a certain parent type        
%       for pp=1:parn(pc)
%           if (data(parpos(pp),cs+1) <= lim(pc,1))
%               parinfo(parpos(pp),2)=1;
%           elseif (data(parpos(pp),cs+1) <= lim(pc,2))
%               parinfo(parpos(pp),2)=2;                
%           elseif (data(parpos(pp),cs+1) <= lim(pc,3))
%               parinfo(parpos(pp),2)=3;                
%           else
%               parinfo(parpos(pp),2)=4;                 
%           end
%       end
%   end

% based on three volume classes, make volume limits (3 volume clases)
                
    lim(1,1:2)=prctile(parval(1:parn(1),1),[33 67]);  
    lim(2,1:2)=prctile(parval(1:parn(2),2),[33 67]);  
    lim(3,1:2)=prctile(parval(1:parn(3),3),[33 67]);   
    lim(4,1:2)=prctile(parval(1:parn(4),4),[33 67]);   
    
    lim;

% add volume class in parinfo, column 2: 1 - smallest,...,5 - highest
    
    for pc=1:4
        parpos=find(parinfo(:,1)==pc); %positions of a certain parent type        
        for pp=1:parn(pc)
            if (data(parpos(pp),cs+1) <= lim(pc,1))
                parinfo(parpos(pp),2)=1;
            elseif (data(parpos(pp),cs+1) <= lim(pc,2))
                parinfo(parpos(pp),2)=2;                
            else
                parinfo(parpos(pp),2)=3;                 
            end
        end
    end  
    
% calculate p-statistics

    for pc=1:4
        ppos=find(parinfo(:,1)==pc); %positions of a certain parent class        
        pinfo=parinfo(ppos,:);   
        
% in the evaluation of probabilities, all volume classes are used
% results in matrix prob        
        
        for vc=1:nvc 
            vpos=find(pinfo(:,2)==vc); %positions of a certain volume class  
            pvinfo=pinfo(vpos,:);   
            pvn=size(pvinfo,1); %number of parents in this pos/vol class
          
            if (pvn==0)
                cs;
                pc;
                vc;
            end
            
            for dt=1:3
                n(cs,pc,vc,dt)=size(find(pvinfo(:,3)==dt),1); %number of parents in this pos/vol/dt class
                prob(cs,pc,vc,dt)=n(cs,pc,vc,dt)/pvn;
            end
            
        end

% in the evaluation of x/x-histograms, all volume classes are pooled
% results in matrix hist

        xpos=find(pinfo(:,3)==3); %positions of x/x-divisions
        xn=size(xpos,1);        
        histo(cs,pc,1:xn)=pinfo(xpos,4);  
        
    end

end

% save observed probabilities for comparison with model probabilities in
% RCD_... (note that scale is reversed in observed and modeled arrays)

pxxo=prob(:,:,:,3);
p01o=prob(:,:,:,1);
p10o=prob(:,:,:,2);

% probabilities - general results

x1v=[ones(nvc,1) (1:nvc)']; %vol class numbers = x in regression (+ a first row of ones, needed for regress in model 1)

scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)])
i=0;
for pc=1:4
    for dt=1:3
        i=i+1;
        subplot(4,3,i);
        pcc=num2str(pc);
        dtc=num2str(dt);   
        for vc=1:nvc %all probabilities are plotted
            pplot(vc,:)=prob(:,pc,vc,dt);
        end
        nsum=sum(n(:,pc,:,dt)); %number of parents in this pos/dt class
        nsumc=num2str(sum(nsum));
        plot(pplot);
        xlim([1 3]);
        set(gca,'XTick',[1 2 3 4]);
        title(['P:' pcc ' D:' dtc ' N:' nsumc]);
        y=mean(pplot,2); %mean probability for each vol in this pos/dt class = y in regression
        mod1=regress(y,x1v); %model 1 is the linear model of mean probabilities
        for vc=1:nvc
            mod1res(vc,1)=mod1(1)+mod1(2)*vc; 
        end
        hold on
        plot(mod1res,'--rs');
        xlabel('Volume Class');
        ylabel('Probability');
        m1slo(pc,dt)=mod1(2); %slope of mean regression, save for model 2 below
    end
end


% probabilities - parameter estimation and plotting

scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)])
i=0;
for pc=1:4
    for dt=1:3
        i=i+1;
        pcc=num2str(pc);
        dtc=num2str(dt);   
        for vc=1:nvc %all probabilities are plotted
            pplot(vc,:)=prob(:,pc,vc,dt);
        end
        nsum=sum(n(:,pc,:,dt)); %number of parents in this pos/dt class
        nsumc=num2str(sum(nsum));
        m1sloc=num2str(m1slo(pc,dt)); %mean regression slope as string
        mod2c=['a+' m1sloc '*x']; %string expression of model 2, i.e. linear regression with fixed mean slope
        mod2=fittype(mod2c); %type of general model to be fitted
        x2v=x1v(:,2); %only the vol class numbers
        for cs=1:nstep
            m2res=fit(x2v,pplot(:,cs),mod2,'Startpoint',0); %linear fit with mean slope to each cascade step 
            m2int(cs,pc,dt)=coeffvalues(m2res); %m2int contains the resulting intercept
        end
        subplot(4,4,i);
        plot(m2int(:,pc,dt));
        xlabel('Cascade Step');
        ylabel('Value of a');
        xlim([1 nstep]);
        title(['Int-variation P:' pcc ' D:' dtc ' Slo:' m1sloc]);
        if (dt==2)
            d01nsum=sum(n(:,pc,:,1),3); %number of "0/1-parents" of this pc, all vcs, for each cs 
            d10nsum=sum(n(:,pc,:,2),3); %number of "1/0-parents" of this pc, all vcs, for each cs 
            d01frac=d01nsum./(d01nsum+d10nsum); %fraction of 0/1-parents of all "non-x/x-parents"
            i=i+1;
            subplot(4,4,i);
            plot(d01frac);
            xlim([1 nstep]);
            xlabel('Cascade Step');
            ylabel('Value of a');
            title(['Fraction 01/(01+10) P:' pcc]);
        end
    end
end

% probabilities - simplified model information 

dum1=cst';
cstlog2=log2(dum1);
x1p=[ones(nstep,1) cstlog2];

% Position class 1 (isolated) 
slom(1,3)=m1slo(1,3); %mean slope of pxx
y=m2int(:,1,3); %associated intercept
mod=regress(y,x1p);
c1(1,3)=mod(1);
c2(1,3)=mod(2);

% Position class 3 (enclosed)
slom(3,3)=m1slo(3,3); %mean slope of pxx
y=m2int(:,3,3); %associated intercept
mod=regress(y,x1p);
c1(3,3)=mod(1);
c2(3,3)=mod(2);

% Position classes 2 and 4 (edge) 
dum4=[m1slo(2,3) m1slo(4,3)]; %pool pxx slopes from 2 and 4
slom(2,3)=mean(dum4); %totalt mean slope of pxx
dum5=[m1slo(2,1) m1slo(4,2)]; %pool p01 and p01 slopes from 2 and 4
slom(2,1)=mean(dum5); %totalt mean slope of p01 and p10
dum2=[m2int(:,2,3) m2int(:,4,3)]; %pool pxx intercepts
y=mean(dum2,2); %mean pxx intercept
mod=regress(y,x1p);
c1(2,3)=mod(1);
c2(2,3)=mod(2);
dum3=[m2int(:,2,1);m2int(:,4,2)]; %pool p01 and p10 intercepts
intm21=mean(dum3); %mean p01 and p10 intercept (constant in the model)

% histograms - shape

figure('Position',[1 1 scrsz(3) scrsz(4)])
i=0;
for pc=1:4
    for cs=1:nstep
        i=i+1;
        subplot(4,nstep,i);
        hplotp=find(histo(cs,pc,:)>0);
        xn=size(hplotp,1);
        clear hplot
        hplot(1:xn)=histo(cs,pc,hplotp);
        for j=1:xn
            hplot(xn+j)=1-hplot(j); %in histo only one weight is given, here the complementary one is added to get symmetry
        end
        if (xn>4) %if less than 5 "unique" weights (i.e. 10 with the complementary), nothing is done (OK?)
            bins=round(1+3.3*log10(2*xn)); %number of bins according to recommended formula
            [nn,bc]=hist(hplot,bins); %create histogram
            bar(bc,nn/sum(nn)); %plot using probabilities instead of number of values (bar)
            afit=betafit(hplot); %fit beta distribution
            ao(cs,pc)=afit(1); %save the a parameter
            bp=betapdf(bc,afit(1),afit(1)); %generate the fitted distribution
            hold on
            plot(bc,bp/sum(bp)); %plot the fit (line)
        else
            hist(0);
        end
        pcc=num2str(pc);
        csc=num2str(cs);
        xnc=num2str(xn);
        xlabel('Weight');
        ylabel('Probability');
        title(['P:' pcc ' CS:' csc ' N:' xnc]);
    end
end

% histograms - a-parameter

figure('Position',[1 1 scrsz(3) scrsz(4)])
i=0;
for pc=1:4
    i=i+1;
    subplot(4,2,i);
    plot(ao(:,pc));
    pcc=num2str(pc);
    xlabel('Cascade Step');
    ylabel('Value of a');
    title(['Beta-a variation lin-log2 P:' pcc]);
    i=i+1;
    subplot(4,2,i);
    balog2(:,pc)=log2(ao(:,pc));
    plot(cstlog2,balog2(:,pc));
    xlabel('Cascade Step');
    ylabel('Value of a');
    title(['Beta-a variation log2-log2 P:' pcc]);    
end

% histograms - simplified model information 

alcs=5; %cascade step below which c3 and c4 is used / above which constant ac is used
alim=x1p(alcs,2); %scale below which c3 and c4 is used / above which constant ac is used

% Position class 1 (isolated) 
mod=regress(balog2(1:alcs,1),x1p(1:alcs,:));
c3(1)=mod(1);
c4(1)=mod(2);
ac(1)=mean(balog2(alcs:nstep,1));
    
% Position class 3 (enclosed) 
mod=regress(balog2(1:alcs,3),x1p(1:alcs,:));
c3(3)=mod(1);
c4(3)=mod(2);
ac(3)=mean(balog2(alcs:nstep,3));

% Position classes 2 and 4 (edge) 
dum6=[balog2(1:alcs,2) balog2(1:alcs,4)];
y=mean(dum6,2);
mod=regress(y,x1p(1:alcs,:));
c3(2)=mod(1);
c4(2)=mod(2);
dum7=[balog2(alcs:nstep,2);balog2(alcs:nstep,4)];
ac(2)=mean(dum7);
save B11.mat
%save 'B.mat' slom c1 c2 intm21 c3 c4 alim ac ao pxxo p01o p10o;