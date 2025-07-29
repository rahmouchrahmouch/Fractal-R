    % This is a general program to perform a Rainfall Cascade Disaggregation
    % The program... 
    
    clear
    close all
    
    % load data, assumed to be in textfile with one column and nr rows, where nr is a number 2^n
    
    %file='C.txt'; %file with data to be disaggregated, unit supposed to be NUMBER OF TIPS
     pfile='B11.mat'; %file with model parameters
    %  load B1.mat
    % set the number of "cascade steps" to be evaluated
    % starting from the highest resolution
    nstep=7; %cascade steps
    tres=675; %final time resolution of the disaggregated data in seconds
    vres=0.0625; %volume resolution of the measurement device
    boxs=[128;64;32;16;8;4;2;1];
    nvc=3; %volume classes 
    boxt=boxs*tres;
    for cs=1:nstep
        cst(cs)=mean(boxt(cs:cs+1));
    end
    dum1=cst';
    cstlog2=log2(dum1);
    
    % Fist make matrices pxxm/p01m/p10m and am of model probabilities and a-values
    
    %  inpname=['C:\Users\Sihem Jebari\Desktop\fractal\Arun_files\',pfile];
    inpname=['D:\test fractal\',pfile];
    %  inpname='B.mat'
    load(inpname);
    
    % Position class 1 (isolated)
    
    for cs=1:nstep
        scale=cstlog2(cs);
        int=c1(1,3)+c2(1,3)*scale;
        if (scale>alim)
            am(cs,1)=2^ac(1);
        else
            am(cs,1)=2^(c3(1)+c4(1)*scale);
        end
        for vc=1:nvc
            pxxm(cs,1,vc)=int+slom(1,3)*vc;
            if (pxxm(cs,1,vc)>1)
                pxxm(cs,1,vc)=1;
            end
            p01m(cs,1,vc)=(1-pxxm(cs,1,vc))/2;
            p10m(cs,1,vc)=(1-pxxm(cs,1,vc))/2;
        end
    end
    
    % Position class 2 (starting)
    
    for cs=1:nstep
        scale=cstlog2(cs);    
        int=c1(2,3)+c2(2,3)*scale;
        if (scale>alim)
            am(cs,2)=2^ac(2);
        else
            am(cs,2)=2^(c3(2)+c4(2)*scale);
        end
        for vc=1:nvc
            pxxm(cs,2,vc)=int+slom(2,3)*vc;
            if (pxxm(cs,2,vc)>1)
                pxxm(cs,2,vc)=1;
            end
            p01m(cs,2,vc)=intm21+slom(2,1)*vc;
            if (p01m(cs,2,vc)>1)
                p01m(cs,2,vc)=1;
            end
            if ((pxxm(cs,2,vc)+p01m(cs,2,vc))>1)
                p01m(cs,2,vc)=1-pxxm(cs,2,vc);
                p10m(cs,2,vc)=0;
            else
                p10m(cs,2,vc)=1-(pxxm(cs,2,vc)+p01m(cs,2,vc));
            end
        end
    end
    
    % Position class 3 (enclosed)
    
    for cs=1:nstep
        scale=cstlog2(cs);      
        int=c1(3,3)+c2(3,3)*scale;
        if (scale>alim)
            am(cs,3)=2^ac(3);
        else
            am(cs,3)=2^(c3(3)+c4(3)*scale);
        end    
        for vc=1:nvc
            pxxm(cs,3,vc)=int+slom(3,3)*vc;
            if (pxxm(cs,3,vc)>1)
                pxxm(cs,3,vc)=1;
            end
            p01m(cs,3,vc)=(1-pxxm(cs,3,vc))/2;
            p10m(cs,3,vc)=(1-pxxm(cs,3,vc))/2;
        end
    end
    
    % Position class 4 (stopping)
    
    for cs=1:nstep
        scale=cstlog2(cs);       
        int=c1(2,3)+c2(2,3)*scale;
        if (scale>alim)
            am(cs,4)=2^ac(2);
        else
            am(cs,4)=2^(c3(2)+c4(2)*scale);
        end    
        for vc=1:nvc
            pxxm(cs,4,vc)=int+slom(2,3)*vc;
            if (pxxm(cs,4,vc)>1)
                pxxm(cs,4,vc)=1;
            end
            p10m(cs,4,vc)=intm21+slom(2,1)*vc;
            if (p10m(cs,4,vc)>1)
                p10m(cs,4,vc)=1;
            end
            if ((pxxm(cs,4,vc)+p10m(cs,4,vc))>1)
                p10m(cs,4,vc)=1-pxxm(cs,4,vc); 
                p01m(cs,4,vc)=0;
            else
                p01m(cs,4,vc)=1-(pxxm(cs,4,vc)+p10m(cs,4,vc));
            end
        end
    end
    
    if pxxm+p01m+p10m ~= 1
        error('sum is not 1')
    end
    %check that the sum of probabilities are always 1
    
    % OBS: when comparing model probabilities (pxxm etc.) with observed
    % probabilites (pxxo etc.), note that the scale is reversed in the
    % different arrays (pxxo start from the smallest scale and pxxm from the
    % largest)
    
    % indname=['C:\Arun Rana\My Dropbox\PhD\Matlab Files and Data\Disaggregation of Rainfall\Mumbai\',file];
    % load Mumbai_24h_in.txt  %file with data to be disaggregated
    % data = Mumbai_24h_in;
    % [nr,nc] = size(data); % nc - columns (1), nr - rows
    
    % indname=['C:\Users\Sihem Jebari\Desktop\sihem\C.dat'];
    load C11.mat
    data=C11;
    [nr,nc] = size(data); % nc - columns (1), nr - rows
    
    
    for nt=1:10
    
    % Disaggregation model
    
    
    for cs=1:nstep
        cs;
    
    % start by creating matrices with "parent information"
    
        npar=nr*(2^(cs-1)); %number of "parents"
        nkid=nr*(2^cs); %number of "kids"
        parn(1:4)=0; %types of parents
        clear parinfo;
        for ap=1:npar
            if data(ap,cs) <0
                data(ap,cs) = 0;
            end
            if (ap==1) %start of series, assume 0 before
                lef=0;
            else
                lef=data(ap-1,cs);
            end
            if (ap==npar) %end of series, assume 0 after
                rig=0;
            else
                rig=data(ap+1,cs);
            end 
            if (data(ap,cs)>0)        
                if ((lef==0) && (rig==0))
                    parn(1)=parn(1)+1;
                    parval(parn(1),1)=data(ap,cs); %matrix with "parent values"
                    parinfo(ap,1)=1; %column 1 "parent types": 1 - isolated, 2 - starting, 3 - enclosed, 4 - stopping
                elseif ((lef==0) && (rig>0))
                    parn(2)=parn(2)+1; 
                    parval(parn(2),2)=data(ap,cs);    
                    parinfo(ap,1)=2;                 
                elseif ((lef>0) && (rig>0))
                    parn(3)=parn(3)+1;
                    parval(parn(3),3)=data(ap,cs);    
                    parinfo(ap,1)=3;                    
                else
                    parn(4)=parn(4)+1;
                    parval(parn(4),4)=data(ap,cs);    
                    parinfo(ap,1)=4;                    
                end
            end
        end
    
        parn;
    
    % based on three volume classes, make volume limits
                    
        vlim(1,1:2)=prctile(parval(1:parn(1),1),[33 67]);  
        vlim(2,1:2)=prctile(parval(1:parn(2),2),[33 67]);  
        vlim(3,1:2)=prctile(parval(1:parn(3),3),[33 67]);   
        vlim(4,1:2)=prctile(parval(1:parn(4),4),[33 67]);   
        
        vlim;
    
    % add volume class in parinfo, column 2: 1 - smallest,...,3 - highest
        
        for pc=1:4
            parpos=find(parinfo(:,1)==pc); %positions of a certain parent type        
            for pp=1:parn(pc)
                if (data(parpos(pp),cs) <= vlim(pc,1))
                    parinfo(parpos(pp),2)=1;
                elseif (data(parpos(pp),cs) <= vlim(pc,2))
                    parinfo(parpos(pp),2)=2;                
                else
                    parinfo(parpos(pp),2)=3;                 
                end
            end
        end  
    
    % Disaggregation start    
        
        for ap=1:npar
    %       cs,ap
            pos1=(2*ap)-1;
            pos2=2*ap;
            if data(ap,cs)==0
                data(pos1,cs+1)=0;
                data(pos2,cs+1)=0;
            else
                pc=parinfo(ap,1); %pc can be 1,2,3, or 4 and not 0
                vc=parinfo(ap,2); %vc can be 1,2,or 3 and not 0
                if (data(ap,cs)==1) %if data is one tip, x/x-div is not possible, therefore rescale p01 and p10 to sum one          
                    p01ms=p01m(cs,pc,vc)/(p01m(cs,pc,vc)+p10m(cs,pc,vc));
                    p10ms=p10m(cs,pc,vc)/(p01m(cs,pc,vc)+p10m(cs,pc,vc));
                    ran=rand;        
                    if (ran<p01ms)
                        data(pos1,cs+1)=0;
                        data(pos2,cs+1)=1;
                    else
                        data(pos1,cs+1)=1;
                        data(pos2,cs+1)=0;
                    end       
                else
                    ran=rand;
                    if (ran<p01m(cs,pc,vc))
                        data(pos1,cs+1)=0;
                        data(pos2,cs+1)=data(ap,cs);
                    elseif (ran<(p01m(cs,pc,vc)+p10m(cs,pc,vc)))
                        data(pos1,cs+1)=data(ap,cs);
                        data(pos2,cs+1)=0;
                    else
                        ran2(1)=gamrnd(am(cs,pc),1);
                        ran2(2)=gamrnd(am(cs,pc),1);
                        ran3=ran2/sum(ran2);
                        data(pos1,cs+1)=round(ran3(1)*data(ap,cs));
                        data(pos2,cs+1)=round(ran3(2)*data(ap,cs));
                        if (data(pos1,cs+1)==0) % to avoid that an x/x-div results in one zero-values
                            data(pos1,cs+1)=1;
                            data(pos2,cs+1)=data(pos2,cs+1)-1;
                        elseif (data(pos2,cs+1)==0) 
                            data(pos2,cs+1)=1;
                            data(pos1,cs+1)=data(pos1,cs+1)-1;
                        end   
                    end
                end
            end
        end    
    end
              
    % Postprocessing to desired time resolution (check that OK with input series)
    % in this case four 675-s values are converted to three 900-s values
    
    j=round(nkid/4);
    sumdiff=0;
    for i=1:j
        startin=(i*4)-3;
        startut=(i*3)-2;
        in(1)=data(startin,nstep+1);
        in(2)=data(startin+1,nstep+1);
        in(3)=data(startin+2,nstep+1);
        in(4)=data(startin+3,nstep+1);
        ut(1)=round(in(1)+(1/3)*in(2));
        ut(2)=round((2/3)*in(2)+(2/3)*in(3));
        ut(3)=round((1/3)*in(3)+in(4));
        diff=sum(in)-sum(ut); %some minor volume may be added or lost in the rounding, this is to keep track of this
        sumdiff=sumdiff+diff;
        out(startut,nt)=ut(1)*vres;
        out(startut+1,nt)=ut(2)*vres;
        out(startut+2,nt)=ut(3)*vres;
    end
    
    sumdiff;
    
    end
    save 'D11.mat' %out -ascii;
    %save 'D' out -ascii;