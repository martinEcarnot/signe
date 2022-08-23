
set_name={'all','parcG','cepG','parcGcepG','parcGcepC','parcGcepS'};

opts = asca('options');
opts.permtest='off';
opts.permfacts='';

ev_all={};
fact_all={};
for i=1:length(set_name)
    i
    load(['SIGNE_asca_',set_name{i},'.mat']);
%     ired=50
%     r=asca(x(1:ired:end,:),fact(1:ired:end,:),opts);
    r=asca(x,fact,opts);
    expl_var=r.XRes.EffectExplVar;
    for j=1:length(name_fact)
        eval(['new_ev=r.X',char(64+j),'.EffectExplVar;'])
        expl_var=[expl_var;new_ev];
    end
ev_all=[ev_all;expl_var];
fact_all=[fact_all;{name_fact}];
end
%save('SIGNE_asca_ev.mat','set_name','ev_all','fact_all');