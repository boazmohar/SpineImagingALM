figure(1), clf
ht = uicontrol('style','text');
fun = ['p=get(1,''currentpoint''); p=p(1,1:2); ' ...
    'set(ht,''pos'',[p 60 20],''string'',num2str(p))'];
    set(1,'buttondownfcn',@(u,evnt)fn_buttonmotion(fun))
