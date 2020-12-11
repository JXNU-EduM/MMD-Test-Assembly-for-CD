function [r1,r2,r3,lr1,lr2,lr3,tt]=CDI(Dis,IB,J,JM)
    tt=zeros(6,1);
    [Inum,K]=size(IB);
    itr=3000;
    nout1=zeros(J,K);
    nout2=zeros(1,J);
    aout1=zeros(J,K);
    aout2=zeros(1,J);
    iout1=zeros(J,K);
    iout2=zeros(1,J);
    tic
    Dj=Dis;
    for i=1:J
        [~,mindex]=max(Dj);
        nout2(i)=mindex;
        nout1(i,:)=IB(mindex,:);
        Dj(mindex)=-1;
    end
    tt(1,1)=toc;
    tt(4,1)=toc;
    tic
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
    beq=[JM(2,:) sum(JM(2,:))];
    lb=zeros(1,Inum);
    ub=ones(1,Inum);
    intcon=(1:Inum);
    f=-Dis;
    options = optimoptions('intlinprog','MaxNodes',itr,'Display','off');
    [x,fval,exitflag,output]= intlinprog(f,intcon,[],[],Aeq,beq,lb,ub,options);%,options
    j=0;
    Iindex=0;
    for i=1:Inum
        Iindex=Iindex+1;
        if x(i)~=0
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
    Aeq=ones(1,Inum);
    beq=J;
    lb=zeros(1,Inum);
    ub=ones(1,Inum);
    intcon=(1:Inum);
    f=-Dis;
    A=-IB';
    options = optimoptions('intlinprog','MaxNodes',itr,'Display','off');
    [x,fval,exitflag,output]= intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);
    j=0;
    Iindex=0;
    for i=1:Inum
        Iindex=Iindex+1;
        if x(i)~=0
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