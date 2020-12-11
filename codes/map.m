function [ amr,pmr ] = map( R,P,prior,Qs,Qn,QnI)
    l=log(P)*R'+log(1-P)*(1-R)';
    C = bsxfun(@times,exp(l),prior);
    [M,I] = max(C,[],1);
    pmr=mean(QnI==I');
    amr=mean(mean(Qn==Qs(I,:)));
end


