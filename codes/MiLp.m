function [r1,r2,r3,lr1,lr2,lr3,tt]=MiLp(Dis,IB,J,JM)
    tt=zeros(6,1);
    itr=3000;
    [Inum,K]=size(IB);
    M=Inum;
    nout1=zeros(J,K);
    nout2=zeros(1,J);
    aout1=zeros(J,K);
    aout2=zeros(1,J);
    iout1=zeros(J,K);
    iout2=zeros(1,J);
    wf2=1;
    tic
    VERD=Dis;
    f=[-mean(VERD,1) wf2*mean(mean(VERD,1))*J];
    intcon=(1:Inum);
    A=[-VERD -ones(size(VERD,1),1)];
    b3=-mean(VERD,2)'*J;
    Aeq=ones(1,Inum+1);
    Aeq(1,Inum+1)=0;
    beq=J;     
    lb=zeros(1,Inum+1);
    lb(1,Inum+1)=-J;
    ub=ones(1,Inum+1);
    ub(1,Inum+1)=J;
    flag=0;
    numofsub=0;
    while flag==0   
        options = optimoptions('intlinprog','MaxNodes',itr,'Display','off');%
        [x1,fval,exitflag,output]= intlinprog(f,intcon,A,b3,Aeq,beq,lb,ub,options);%,options
        if exitflag<1 
            b3=b3*0.9;
            numofsub=numofsub+1;
        else
            flag=1;
            break;
        end
    end
     j=0;
     Iindex=0;
     for i=1:Inum
         Iindex=Iindex+1;
         if x1(i)~=0
             j=j+1;
             nout1(j,:)=IB(i,:);
             nout2(j)=Iindex;
         end
     end
     tt(1,1)=toc;
     tt(4,1)=toc;
     tic
     VERD=Dis;
     f=[-mean(VERD,1) wf2*mean(mean(VERD,1))*J];
     A=[-VERD -ones(size(VERD,1),1)];
     b3=-mean(VERD,2)'*J;
     SumAttr=sum(IB,2);
     Aeq=zeros(3,Inum);
     Aeq(4,:)=ones(1,Inum);
     for i=1:Inum
         if SumAttr(i)==JM(1,1)
             Aeq(1,i)=1;
         elseif SumAttr(i)==JM(1,2)
             Aeq(2,i)=1;
         elseif SumAttr(i)==JM(1,3)
             Aeq(3,i)=1; 
         end
     end
     Aeq(:,Inum+1)=0;
     beq=[JM(2,:) J];
     lb(1,Inum+1)=0;
     ub(1,Inum+1)=J;  
     flag=0;
     numofsub=0;
     while flag==0 
         options = optimoptions('intlinprog','MaxNodes',itr,'Display','off');%
         [x2,fval,exitflag,output]= intlinprog(f,intcon,A,b3,Aeq,beq,lb,ub,options);    
         if exitflag<1 
             b3=b3*0.9;
             numofsub=numofsub+1;
         else
             flag=1;
             break;
         end
     end
     j=0;
     Iindex=0;
     for i=1:Inum
         Iindex=Iindex+1;
         if x2(i)~=0
             j=j+1;
             iout1(j,:)=IB(i,:);
             iout2(j)=Iindex;
         end
     end
    tt(2,1)=toc;
    tt(5,1)=toc;
    tic
    if K==4
       b=-7*ones(1,K);
    elseif K==8
        b=-3*ones(1,K);
    end
    VERD=Dis;
    f=[-mean(VERD,1) wf2*mean(mean(VERD,1))*J];
    A=[-VERD -ones(size(VERD,1),1)];
    b3=-mean(VERD,2)'*J;
    IBMatr=IB;
    IBMatr(Inum+1,:)=0;
    IBMatr=-IBMatr';
    A=[A;IBMatr];
    b2=[b3,b];
    Aeq=ones(1,Inum+1);
    Aeq(1,Inum+1)=0;
    beq=J; 
    lb(1,Inum+1)=-J;
    ub(1,Inum+1)=J;
    flag=0;
    numofsub=0;
    while flag==0   
        options = optimoptions('intlinprog','MaxNodes',itr,'Display','off');%
        [x3,fval,exitflag,output]= intlinprog(f,intcon,A,b2,Aeq,beq,lb,ub,options);%,options
        if exitflag<1 
            b3=b3*0.95;
            b2=[b3,b,maxb];
            numofsub=numofsub+1;
        else
            flag=1;
            break;
        end
    end
    Iindex=0;
    j=0;
    for i=1:Inum
        Iindex=Iindex+1;
        if x3(i)~=0
            j=j+1;
            aout1(j,:)=IB(i,:);
            aout2(j)=Iindex;
        end
    end
    tt(3,1)=toc;
    tt(6,1)=toc;
    r1=nout1;
    lr1=nout2;
    r2=iout1;
    lr2=iout2;
    r3=aout1;
    lr3=aout2;
end