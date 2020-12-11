function result=ModelParameter(model,An,In)
    if strcmp('DINA',model)||strcmp('DINO',model)
        l=0.05;h=0.40;
        sj=l+(h-l)*rand(1,In);
        gj=l+(h-l)*rand(1,In);
         result=[sj;gj];
    elseif strcmp('R-RUM',model)
        pil=0.75;pih=0.95;
        rl=0.2;rh=0.95;
        pi=pil+(pih-pil)*rand(1,In);
        r=rl+(rh-rl)*rand(An,In);
        result=[pi;r];
    else
        result=-1;
    end
end