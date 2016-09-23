function [H, HI, D1, D2, D4, e_l, e_r, M4, d2_l, d2_r, d3_l, d3_r, d1_l, d1_r] = d4_variable_6_2(m,h)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 6:te ordn. SBP Finita differens         %%%
    %%% operatorer med diagonal norm            %%%
    %%% Extension to variable koeff             %%%
    %%%                                         %%%
    %%% H           (Normen)                    %%%
    %%% D1=H^(-1)Q  (approx f?rsta derivatan)   %%%
    %%% D2          (approx andra derivatan)    %%%
    %%% D2=HI*(R+C*D*S                          %%%
    %%%                                         %%%
    %%% R=-D1'*H*C*D1-RR                        %%%
    %%%                                         %%%
    %%% RR ?r dissipation)                      %%%
    %%% Dissipationen uppbyggd av D4:           %%%
    %%% DI=D4*B*H*D4                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % H?r med 6 RP ist?llet f?r 8 f?r D4 operatorn, dock samma randderivator
    % Denna ?r noggrannare, och har 2a ordningens randdslutning och b?r ge 6te
    % ordningens konvergens. Hade dock ingen fri parameter att optimera


    H = diag(ones(m,1),0);
    H(1:6,1:6) = [
        0.181e3/0.576e3 0 0 0 0 0;
        0 0.1343e4/0.960e3 0 0 0 0;
        0 0 0.293e3/0.480e3 0 0 0;
        0 0 0 0.1811e4/0.1440e4 0 0;
        0 0 0 0 0.289e3/0.320e3 0;
        0 0 0 0 0 0.65e2/0.64e2;
    ];
    H(m-5:m,m-5:m) = fliplr(flipud(H(1:6,1:6)));


    e_1 = zeros(m,1);e_1(1) = 1;
    e_m = zeros(m,1);e_m(m) = 1;

    S_U = [-0.137e3/0.60e2 5 -5 0.10e2/0.3e1 -0.5e1/0.4e1 0.1e1/0.5e1;]/h;
    S_1 = zeros(1,m);
    S_1(1:6) = S_U;
    S_m = zeros(1,m);
    S_m(m-5:m) = fliplr(-S_U);



    %DS = zeros(m,m);
    %DS(1,1:5) = -[-25/12, 4, -3, 4/3, -1/4];
    %DS(m,m-4:m) = fliplr(-[-25/12, 4, -3, 4/3, -1/4]);
    %DS = diag(c)*DS/h;


    H = h*H;
    HI = inv(H);


    % M = zeros(m,m);
    %
    % for i=6:m-5
    % M(i,i-3:i+3)=[c(i-2)/0.40e2 + c(i-1)/0.40e2 - 0.11e2/0.360e3 * c(i-3) - 0.11e2/0.360e3 * c(i) c(i-3)/0.20e2 - 0.3e1/0.10e2 * c(i-1) + c(i+1)/0.20e2 + 0.7e1/0.40e2 * c(i) + 0.7e1/0.40e2 * c(i-2) -c(i-3)/0.40e2 - 0.3e1/0.10e2 * c(i-2) - 0.3e1/0.10e2 * c(i+1) - c(i+2)/0.40e2 - 0.17e2/0.40e2 * c(i) - 0.17e2/0.40e2 * c(i-1) c(i-3)/0.180e3 + c(i-2)/0.8e1 + 0.19e2/0.20e2 * c(i-1) + 0.19e2/0.20e2 * c(i+1) + c(i+2)/0.8e1 + c(i+3)/0.180e3 + 0.101e3/0.180e3 * c(i) -c(i-2)/0.40e2 - 0.3e1/0.10e2 * c(i-1) - 0.3e1/0.10e2 * c(i+2) - c(i+3)/0.40e2 - 0.17e2/0.40e2 * c(i) - 0.17e2/0.40e2 * c(i+1) c(i-1)/0.20e2 - 0.3e1/0.10e2 * c(i+1) + c(i+3)/0.20e2 + 0.7e1/0.40e2 * c(i) + 0.7e1/0.40e2 * c(i+2) c(i+1)/0.40e2 + c(i+2)/0.40e2 - 0.11e2/0.360e3 * c(i) - 0.11e2/0.360e3 * c(i+3);];
    % end
    %
    %
    % M(1:9,1:9)=[
    %     0.7912667594695582093926295e0 * c(1) + 0.2968472090638000742888467e0 * c(2) + 0.3185519088796429015220016e-2 * c(3) + 0.1632404042590951953384672e-1 * c(4) + 0.3160302244094415087693968e-1 * c(5) + 0.3167964748016105299646518e-1 * c(6) + 0.3148577733947253920469418e-1 * c(7) -0.1016689339350338144430605e1 * c(1) - 0.2845627370491611369031341e-1 * c(3) - 0.4128029838349298819825156e-1 * c(4) - 0.1392281451620140507549866e0 * c(5) - 0.1195777325611201766551392e0 * c(6) - 0.1194267756529333410855186e0 * c(7) 0.7075642937243715046279337e-1 * c(1) - 0.1845476106024151050283847e0 * c(2) - 0.4364163147111892346990101e-1 * c(4) + 0.2432367907207732460874765e0 * c(5) + 0.1582127073537215443965653e0 * c(6) + 0.1602348578364786307613271e0 * c(7) 0.2251991532891353212689574e0 * c(1) - 0.1662748711097054895317080e0 * c(2) + 0.2710530961648671297733465e-1 * c(3) - 0.1916646185968439909125616e0 * c(5) - 0.7684117160199014594442072e-1 * c(6) - 0.8219586949831697575883635e-1 * c(7) -0.5224403464202056316702078e-1 * c(1) + 0.4440063948509876221050939e-1 * c(2) - 0.1023976547309387874453988e-2 * c(3) + 0.7403484645316174090533193e-1 * c(4) + 0.1241625568998496895352046e-1 * c(6) + 0.7188652847892601282652228e-1 * c(5) + 0.1379362997104735503447960e-1 * c(7) -0.1828896813877197352675410e-1 * c(1) + 0.9574633163221758060736592e-2 * c(2) - 0.8105784530576404277872603e-3 * c(3) - 0.7348845587775519698437916e-2 * c(4) + 0.1063601949723906997026904e-1 * c(5) - 0.1315967038382618382356495e-1 * c(6) - 0.2117936478838753524581943e-1 * c(7) 0.1911888563316170927411831e-2 * c(4) - 0.4068130355529149936100229e-1 * c(5) + 0.1319674981073749167009902e-1 * c(6) + 0.2557266518123783676349144e-1 * c(7) 0.1559652871136785763960685e-1 * c(5) - 0.6486184157331537899459796e-2 * c(6) - 0.9110344554036319740147054e-2 * c(7) 0.5593983696629863059347067e-3 * c(6) - 0.1384822535100796372263822e-2 * c(5) + 0.8254241654378100663291154e-3 * c(7);
    %     -0.1016689339350338144430605e1 * c(1) - 0.2845627370491611369031341e-1 * c(3) - 0.4128029838349298819825156e-1 * c(4) - 0.1392281451620140507549866e0 * c(5) - 0.1195777325611201766551392e0 * c(6) - 0.1194267756529333410855186e0 * c(7) 0.1306332157111667628555907e1 * c(1) + 0.2542001760457345743492403e0 * c(3) + 0.1043897828092562609502636e0 * c(4) + 0.6672328021032112950919876e0 * c(5) + 0.4681819359722749441073885e0 * c(6) + 0.4676415410195836920069412e0 * c(7) -0.9091410269992464604926176e-1 * c(1) + 0.1103611313171476425250639e0 * c(4) - 0.1290397544997518887000350e1 * c(5) - 0.6639605248735044787146222e0 * c(6) - 0.6615974464005206184151509e0 * c(7) -0.2893557395653431666593814e0 * c(1) - 0.2421320004064592721552708e0 * c(3) + 0.1187670255028031027693374e1 * c(5) + 0.3956598149904136332753521e0 * c(6) + 0.3860048921755800000681479e0 * c(7) 0.6712774475803763988977355e-1 * c(1) + 0.9147192682075630179962131e-2 * c(3) - 0.1872196143003808021730728e0 * c(4) - 0.1319358558853174530078498e0 * c(6) - 0.4871575736811911887376923e0 * c(5) - 0.1047516312275448138054418e0 * c(7) 0.2349927974590068869356781e-1 * c(1) + 0.7240905383565181316381731e-2 * c(3) + 0.1858378996391679448655070e-1 * c(4) - 0.9289616133938676174345208e-1 * c(5) + 0.1223513270418807666970488e0 * c(6) + 0.1113520320436295033894092e0 * c(7) -0.4834791406446907590553793e-2 * c(4) + 0.2310683832687820403062715e0 * c(5) - 0.1080774142196007991746827e0 * c(6) - 0.1181561776427343335410349e0 * c(7) -0.8368141434403455353724691e-1 * c(5) + 0.4093499466767054661591066e-1 * c(6) + 0.4274641967636400692133626e-1 * c(7) -0.3576545132696983143406173e-2 * c(6) + 0.7389399124121078682094445e-2 * c(5) - 0.3812853991424095538688273e-2 * c(7);
    %     0.7075642937243715046279337e-1 * c(1) - 0.1845476106024151050283847e0 * c(2) - 0.4364163147111892346990101e-1 * c(4) + 0.2432367907207732460874765e0 * c(5) + 0.1582127073537215443965653e0 * c(6) + 0.1602348578364786307613271e0 * c(7) -0.9091410269992464604926176e-1 * c(1) + 0.1103611313171476425250639e0 * c(4) - 0.1290397544997518887000350e1 * c(5) - 0.6639605248735044787146222e0 * c(6) - 0.6615974464005206184151509e0 * c(7) 0.6327161147136873807796515e-2 * c(1) + 0.1147318200715868527529827e0 * c(2) + 0.1166740554279680007487795e0 * c(4) + 0.2766610808285444037240703e1 * c(5) + 0.1070920689960817104203947e1 * c(6) + 0.1013161391032973057171717e1 * c(7) 0.2013769413884797246646959e-1 * c(1) + 0.1033717994630886401730470e0 * c(2) - 0.2913221621151742724258117e1 * c(5) - 0.8755807343482262259774782e0 * c(6) - 0.6909957183488812426508351e0 * c(7) -0.4671751091575462868310238e-2 * c(1) - 0.2760353365637712827793337e-1 * c(2) - 0.1979290298620869974478871e0 * c(4) + 0.5402985338373433052255418e0 * c(6) + 0.1239177593031911077924537e1 * c(5) + 0.2628038050247358227280031e0 * c(7) -0.1635430866921887819487473e-2 * c(1) - 0.5952475275883259619711594e-2 * c(2) + 0.1964682777744275219350831e-1 * c(4) + 0.3236640012639046600590714e0 * c(5) - 0.4659516693228870973898560e0 * c(6) - 0.2217272720941736859420432e0 * c(7) -0.5111353189352474549563559e-2 * c(4) - 0.5355878163774754346032096e0 * c(5) + 0.3328335104489738933610597e0 * c(6) + 0.2078656591178540157917135e0 * c(7) 0.1824328174134289562208038e0 * c(5) - 0.1059816030196818445908057e0 * c(6) - 0.7645121439374711162999809e-1 * c(7) 0.9209089963443799485648361e-2 * c(6) - 0.1591502818872493167091475e-1 * c(5) + 0.6705938225281132185266388e-2 * c(7);
    %     0.2251991532891353212689574e0 * c(1) - 0.1662748711097054895317080e0 * c(2) + 0.2710530961648671297733465e-1 * c(3) - 0.1916646185968439909125616e0 * c(5) - 0.7684117160199014594442072e-1 * c(6) - 0.8219586949831697575883635e-1 * c(7) -0.2893557395653431666593814e0 * c(1) - 0.2421320004064592721552708e0 * c(3) + 0.1187670255028031027693374e1 * c(5) + 0.3956598149904136332753521e0 * c(6) + 0.3860048921755800000681479e0 * c(7) 0.2013769413884797246646959e-1 * c(1) + 0.1033717994630886401730470e0 * c(2) - 0.2913221621151742724258117e1 * c(5) - 0.8755807343482262259774782e0 * c(6) - 0.6909957183488812426508351e0 * c(7) 0.6409299775987186986730499e-1 * c(1) + 0.9313657638804699948929701e-1 * c(2) + 0.2306367624634749229113646e0 * c(3) + 0.3689440308283716620260816e1 * c(5) + 0.1190550338687608873798462e1 * c(6) + 0.5912479546888856519443605e0 * c(7) -0.1486895819265604128572498e-1 * c(1) - 0.2487040599390160764166412e-1 * c(2) - 0.8712928907711754187084757e-2 * c(3) - 0.1263507837371824205693950e1 * c(6) - 0.3058317397843997326920898e0 * c(7) - 0.1470691926045802954795783e1 * c(5) -0.5205147429855955657625694e-2 * c(1) - 0.5363098747528542488971874e-2 * c(2) - 0.6897142765790609546343709e-2 * c(3) - 0.7857524521667450101721993e0 * c(5) + 0.2291148005423734600066709e0 * c(7) + 0.9977064356292750529201981e0 * c(6) 0.6697297488067662265210608e0 * c(5) - 0.5013247356072127938999311e0 * c(6) - 0.1795161243106645437322408e0 * c(7) -0.2022909060111751565150958e0 * c(5) + 0.1453421858063658498587377e0 * c(6) + 0.5694872020480930665635812e-1 * c(7) -0.1200429618441003833696998e-1 * c(6) - 0.4776915669385923841535432e-2 * c(7) + 0.1678121185379596217850541e-1 * c(5);
    %     -0.5224403464202056316702078e-1 * c(1) + 0.4440063948509876221050939e-1 * c(2) - 0.1023976547309387874453988e-2 * c(3) + 0.7403484645316174090533193e-1 * c(4) + 0.1241625568998496895352046e-1 * c(6) + 0.7188652847892601282652228e-1 * c(5) + 0.1379362997104735503447960e-1 * c(7) 0.6712774475803763988977355e-1 * c(1) + 0.9147192682075630179962131e-2 * c(3) - 0.1872196143003808021730728e0 * c(4) - 0.1319358558853174530078498e0 * c(6) - 0.4871575736811911887376923e0 * c(5) - 0.1047516312275448138054418e0 * c(7) -0.4671751091575462868310238e-2 * c(1) - 0.2760353365637712827793337e-1 * c(2) - 0.1979290298620869974478871e0 * c(4) + 0.5402985338373433052255418e0 * c(6) + 0.1239177593031911077924537e1 * c(5) + 0.2628038050247358227280031e0 * c(7) -0.1486895819265604128572498e-1 * c(1) - 0.2487040599390160764166412e-1 * c(2) - 0.8712928907711754187084757e-2 * c(3) - 0.1263507837371824205693950e1 * c(6) - 0.3058317397843997326920898e0 * c(7) - 0.1470691926045802954795783e1 * c(5) 0.3449455095910233625229891e-2 * c(1) + 0.6641183499427826101618457e-2 * c(2) + 0.3291545083271862858501887e-3 * c(3) + 0.3357721707576477199985656e0 * c(4) + 0.2096413329579026439044119e1 * c(6) + 0.2317323204183126854954203e0 * c(7) + 0.6107825764368264576481962e-2 * c(8) + 0.7109125850683376695640722e0 * c(5) 0.1207544072304193806052558e-2 * c(1) + 0.1432116665752147607469646e-2 * c(2) + 0.2605582646183255957264249e-3 * c(3) - 0.3332941113251635390801278e-1 * c(4) - 0.2808241697385532683612407e0 * c(7) - 0.2720908083525083608370563e-1 * c(8) + 0.1045865435682921987447929e0 * c(5) - 0.1348436986667115543203552e1 * c(6) 0.8671038084174692625075159e-2 * c(4) + 0.1736073411355428563685818e0 * c(6) + 0.5331362125287625412555844e-1 * c(8) - 0.2424935262404526301801157e0 * c(5) + 0.1569015257678588270609004e0 * c(7) -0.8631683980217122275970376e-1 * c(6) + 0.2698842360470999243492629e-1 * c(7) + 0.8098194147715651085292754e-1 * c(5) - 0.3276463639080639163926118e-1 * c(8) 0.7462059484530855073291365e-2 * c(6) - 0.8121640361668678949573496e-3 * c(7) + 0.5522702088127090209264064e-3 * c(8) - 0.7202165657176696199260422e-2 * c(5);
    %     -0.1828896813877197352675410e-1 * c(1) + 0.9574633163221758060736592e-2 * c(2) - 0.8105784530576404277872603e-3 * c(3) - 0.7348845587775519698437916e-2 * c(4) + 0.1063601949723906997026904e-1 * c(5) - 0.1315967038382618382356495e-1 * c(6) - 0.2117936478838753524581943e-1 * c(7) 0.2349927974590068869356781e-1 * c(1) + 0.7240905383565181316381731e-2 * c(3) + 0.1858378996391679448655070e-1 * c(4) - 0.9289616133938676174345208e-1 * c(5) + 0.1223513270418807666970488e0 * c(6) + 0.1113520320436295033894092e0 * c(7) -0.1635430866921887819487473e-2 * c(1) - 0.5952475275883259619711594e-2 * c(2) + 0.1964682777744275219350831e-1 * c(4) + 0.3236640012639046600590714e0 * c(5) - 0.4659516693228870973898560e0 * c(6) - 0.2217272720941736859420432e0 * c(7) -0.5205147429855955657625694e-2 * c(1) - 0.5363098747528542488971874e-2 * c(2) - 0.6897142765790609546343709e-2 * c(3) - 0.7857524521667450101721993e0 * c(5) + 0.2291148005423734600066709e0 * c(7) + 0.9977064356292750529201981e0 * c(6) 0.1207544072304193806052558e-2 * c(1) + 0.1432116665752147607469646e-2 * c(2) + 0.2605582646183255957264249e-3 * c(3) - 0.3332941113251635390801278e-1 * c(4) - 0.2808241697385532683612407e0 * c(7) - 0.2720908083525083608370563e-1 * c(8) + 0.1045865435682921987447929e0 * c(5) - 0.1348436986667115543203552e1 * c(6) 0.4227226173449345042468960e-3 * c(1) + 0.3088241944378964404772302e-3 * c(2) + 0.2062575706647430620228133e-3 * c(3) + 0.3308343404200968256656458e-2 * c(4) + 0.5828047016405001815804837e0 * c(5) + 0.8054174220366215473556835e0 * c(7) + 0.1338363233410033443348225e0 * c(8) + 0.5555555555555555555555556e-2 * c(9) + 0.1190362071861893051132274e1 * c(6) -0.8607044252686413302647675e-3 * c(4) - 0.1748074708673904989293256e0 * c(5) - 0.3132544850115050165022338e0 * c(8) - 0.2500000000000000000000000e-1 * c(9) - 0.3169166305310429271303167e0 * c(7) - 0.6691607091647929161078591e0 * c(6) 0.3354661791693352108660900e-1 * c(5) - 0.3343620022386971405018586e0 * c(7) + 0.5000000000000000000000000e-1 * c(9) + 0.2169790609807602750804271e0 * c(6) + 0.1838363233410033443348225e0 * c(8) 0.2912518476823004642951502e-1 * c(7) + 0.2279091916474916391629437e-1 * c(8) - 0.3068985997518740530511593e-1 * c(6) - 0.1781799513347360596249022e-2 * c(5) - 0.3055555555555555555555556e-1 * c(9);
    %     0.1911888563316170927411831e-2 * c(4) - 0.4068130355529149936100229e-1 * c(5) + 0.1319674981073749167009902e-1 * c(6) + 0.2557266518123783676349144e-1 * c(7) -0.4834791406446907590553793e-2 * c(4) + 0.2310683832687820403062715e0 * c(5) - 0.1080774142196007991746827e0 * c(6) - 0.1181561776427343335410349e0 * c(7) -0.5111353189352474549563559e-2 * c(4) - 0.5355878163774754346032096e0 * c(5) + 0.3328335104489738933610597e0 * c(6) + 0.2078656591178540157917135e0 * c(7) 0.6697297488067662265210608e0 * c(5) - 0.5013247356072127938999311e0 * c(6) - 0.1795161243106645437322408e0 * c(7) 0.8671038084174692625075159e-2 * c(4) + 0.1736073411355428563685818e0 * c(6) + 0.5331362125287625412555844e-1 * c(8) - 0.2424935262404526301801157e0 * c(5) + 0.1569015257678588270609004e0 * c(7) -0.8607044252686413302647675e-3 * c(4) - 0.1748074708673904989293256e0 * c(5) - 0.3132544850115050165022338e0 * c(8) - 0.2500000000000000000000000e-1 * c(9) - 0.3169166305310429271303167e0 * c(7) - 0.6691607091647929161078591e0 * c(6) 0.2239223735771599178951297e-3 * c(4) + 0.1275437785430956673825710e0 * c(5) + 0.1011699483929608164601067e1 * c(6) + 0.9698817275172575247533506e0 * c(8) + 0.1250000000000000000000000e0 * c(9) + 0.5555555555555555555555556e-2 * c(10) + 0.4823177543031281500117826e0 * c(7) -0.3784113973033012949863031e-1 * c(5) - 0.2997556885134827361576001e0 * c(6) - 0.3000000000000000000000000e0 * c(9) - 0.2500000000000000000000000e-1 * c(10) - 0.3991486867446821178415359e0 * c(7) - 0.4382544850115050165022338e0 * c(8) 0.4698146218022683933926520e-1 * c(6) - 0.2966863787471237458744416e0 * c(8) + 0.5000000000000000000000000e-1 * c(10) + 0.1716355704146006481727960e0 * c(7) + 0.3069346152296258362380356e-2 * c(5) + 0.1750000000000000000000000e0 * c(9);
    %     0.1559652871136785763960685e-1 * c(5) - 0.6486184157331537899459796e-2 * c(6) - 0.9110344554036319740147054e-2 * c(7) -0.8368141434403455353724691e-1 * c(5) + 0.4093499466767054661591066e-1 * c(6) + 0.4274641967636400692133626e-1 * c(7) 0.1824328174134289562208038e0 * c(5) - 0.1059816030196818445908057e0 * c(6) - 0.7645121439374711162999809e-1 * c(7) -0.2022909060111751565150958e0 * c(5) + 0.1453421858063658498587377e0 * c(6) + 0.5694872020480930665635812e-1 * c(7) -0.8631683980217122275970376e-1 * c(6) + 0.2698842360470999243492629e-1 * c(7) + 0.8098194147715651085292754e-1 * c(5) - 0.3276463639080639163926118e-1 * c(8) 0.3354661791693352108660900e-1 * c(5) - 0.3343620022386971405018586e0 * c(7) + 0.5000000000000000000000000e-1 * c(9) + 0.2169790609807602750804271e0 * c(6) + 0.1838363233410033443348225e0 * c(8) -0.3784113973033012949863031e-1 * c(5) - 0.2997556885134827361576001e0 * c(6) - 0.3000000000000000000000000e0 * c(9) - 0.2500000000000000000000000e-1 * c(10) - 0.3991486867446821178415359e0 * c(7) - 0.4382544850115050165022338e0 * c(8) 0.1230328942716804455358698e-1 * c(5) + 0.1183647529645898332481833e0 * c(6) + 0.9410511898227943334189628e0 * c(7) + 0.9500000000000000000000000e0 * c(9) + 0.1250000000000000000000000e0 * c(10) + 0.5555555555555555555555556e-2 * c(11) + 0.5699474344521144554459336e0 * c(8) -0.2308067892671916339568942e-1 * c(6) - 0.2986625053775149497180439e0 * c(7) - 0.3000000000000000000000000e0 * c(10) - 0.2500000000000000000000000e-1 * c(11) - 0.1047734860515050802561078e-2 * c(5) - 0.4272090808352508360837056e0 * c(8) - 0.4250000000000000000000000e0 * c(9);
    %     0.5593983696629863059347067e-3 * c(6) - 0.1384822535100796372263822e-2 * c(5) + 0.8254241654378100663291154e-3 * c(7) -0.3576545132696983143406173e-2 * c(6) + 0.7389399124121078682094445e-2 * c(5) - 0.3812853991424095538688273e-2 * c(7) 0.9209089963443799485648361e-2 * c(6) - 0.1591502818872493167091475e-1 * c(5) + 0.6705938225281132185266388e-2 * c(7) -0.1200429618441003833696998e-1 * c(6) - 0.4776915669385923841535432e-2 * c(7) + 0.1678121185379596217850541e-1 * c(5) 0.7462059484530855073291365e-2 * c(6) - 0.8121640361668678949573496e-3 * c(7) + 0.5522702088127090209264064e-3 * c(8) - 0.7202165657176696199260422e-2 * c(5) 0.2912518476823004642951502e-1 * c(7) + 0.2279091916474916391629437e-1 * c(8) - 0.3068985997518740530511593e-1 * c(6) - 0.1781799513347360596249022e-2 * c(5) - 0.3055555555555555555555556e-1 * c(9) 0.4698146218022683933926520e-1 * c(6) - 0.2966863787471237458744416e0 * c(8) + 0.5000000000000000000000000e-1 * c(10) + 0.1716355704146006481727960e0 * c(7) + 0.3069346152296258362380356e-2 * c(5) + 0.1750000000000000000000000e0 * c(9) -0.2308067892671916339568942e-1 * c(6) - 0.2986625053775149497180439e0 * c(7) - 0.3000000000000000000000000e0 * c(10) - 0.2500000000000000000000000e-1 * c(11) - 0.1047734860515050802561078e-2 * c(5) - 0.4272090808352508360837056e0 * c(8) - 0.4250000000000000000000000e0 * c(9) 0.5139370221149109977041877e-2 * c(6) + 0.1247723215009422001393184e0 * c(7) + 0.9505522702088127090209264e0 * c(8) + 0.9500000000000000000000000e0 * c(10) + 0.1250000000000000000000000e0 * c(11) + 0.5555555555555555555555556e-2 * c(12) + 0.9159362465153641826887659e-4 * c(5) + 0.5611111111111111111111111e0 * c(9);
    % ];
    %
    % M(m-8:m,m-8:m)=[
    %     0.5555555555555555555555556e-2 * c(m-11) + 0.1250000000000000000000000e0 * c(m-10) + 0.9500000000000000000000000e0 * c(m-9) + 0.9505522702088127090209264e0 * c(m-7) + 0.1247205076844361998744053e0 * c(m-6) + 0.5139370221149109977041877e-2 * c(m-5) + 0.5611111111111111111111111e0 * c(m-8) + 0.1434074411575366831819799e-3 * c(m-4) -0.2500000000000000000000000e-1 * c(m-10) - 0.3000000000000000000000000e0 * c(m-9) - 0.2980649679116425253322056e0 * c(m-6) - 0.2308067892671916339568942e-1 * c(m-5) - 0.4250000000000000000000000e0 * c(m-8) - 0.4272090808352508360837056e0 * c(m-7) - 0.1645272326387475188399322e-2 * c(m-4) 0.5000000000000000000000000e-1 * c(m-9) - 0.2966863787471237458744416e0 * c(m-7) + 0.4698146218022683933926520e-1 * c(m-5) + 0.1750000000000000000000000e0 * c(m-8) + 0.1700291833903489463825077e0 * c(m-6) + 0.4675733176547960152668626e-2 * c(m-4) 0.2279091916474916391629437e-1 * c(m-7) + 0.3097763128598982561225538e-1 * c(m-6) - 0.3055555555555555555555556e-1 * c(m-8) - 0.3068985997518740530511593e-1 * c(m-5) - 0.3634246031107139778989373e-2 * c(m-4) 0.5522702088127090209264064e-3 * c(m-7) - 0.3265435411305071914756373e-2 * c(m-6) + 0.7462059484530855073291365e-2 * c(m-5) - 0.4748894282038492179461399e-2 * c(m-4) 0.6272075574042975468177820e-3 * c(m-6) - 0.1200429618441003833696998e-1 * c(m-5) + 0.1137708862700574079015220e-1 * c(m-4) 0.9209089963443799485648361e-2 * c(m-5) - 0.3129629392354775191148163e-3 * c(m-6) - 0.8896127024208321966533544e-2 * c(m-4) -0.3576545132696983143406173e-2 * c(m-5) + 0.4335019854436220306755673e-3 * c(m-6) + 0.3143043147253361112730605e-2 * c(m-4) 0.5593983696629863059347067e-3 * c(m-5) - 0.1446656414398166805849327e-3 * c(m-6) - 0.4147327282231696253497740e-3 * c(m-4);
    %     -0.2500000000000000000000000e-1 * c(m-10) - 0.3000000000000000000000000e0 * c(m-9) - 0.2980649679116425253322056e0 * c(m-6) - 0.2308067892671916339568942e-1 * c(m-5) - 0.4250000000000000000000000e0 * c(m-8) - 0.4272090808352508360837056e0 * c(m-7) - 0.1645272326387475188399322e-2 * c(m-4) 0.5555555555555555555555556e-2 * c(m-10) + 0.1250000000000000000000000e0 * c(m-9) + 0.9500000000000000000000000e0 * c(m-8) + 0.9341601509609901526962449e0 * c(m-6) + 0.1183647529645898332481833e0 * c(m-5) + 0.1919432828897222527630486e-1 * c(m-4) + 0.5699474344521144554459336e0 * c(m-7) -0.2500000000000000000000000e-1 * c(m-9) - 0.3000000000000000000000000e0 * c(m-8) - 0.2997556885134827361576001e0 * c(m-5) - 0.5636663150858098975790317e-1 * c(m-4) - 0.4382544850115050165022338e0 * c(m-7) - 0.3806231949664312575822630e0 * c(m-6) 0.5000000000000000000000000e-1 * c(m-8) - 0.3557251496099816106154206e0 * c(m-6) + 0.5490976528821799120017102e-1 * c(m-4) + 0.1838363233410033443348225e0 * c(m-7) + 0.2169790609807602750804271e0 * c(m-5) 0.5528052133944605740009217e-1 * c(m-6) - 0.8631683980217122275970376e-1 * c(m-5) - 0.3276463639080639163926118e-1 * c(m-7) + 0.5268984374242044588776166e-1 * c(m-4) -0.5373770512016897565958305e-2 * c(m-6) + 0.1453421858063658498587377e0 * c(m-5) - 0.1399684152943489522927794e0 * c(m-4) -0.1059816030196818445908057e0 * c(m-5) + 0.1014880675788250237247178e0 * c(m-4) + 0.4493535440856820866087846e-2 * c(m-6) 0.4093499466767054661591066e-1 * c(m-5) - 0.3471075437892810033585296e-1 * c(m-4) - 0.6224240288742446280057699e-2 * c(m-6) -0.6486184157331537899459796e-2 * c(m-5) + 0.4409068609809831485979484e-2 * c(m-4) + 0.2077115547521706413480312e-2 * c(m-6);
    %     0.5000000000000000000000000e-1 * c(m-9) - 0.2966863787471237458744416e0 * c(m-7) + 0.4698146218022683933926520e-1 * c(m-5) + 0.1750000000000000000000000e0 * c(m-8) + 0.1700291833903489463825077e0 * c(m-6) + 0.4675733176547960152668626e-2 * c(m-4) -0.2500000000000000000000000e-1 * c(m-9) - 0.3000000000000000000000000e0 * c(m-8) - 0.2997556885134827361576001e0 * c(m-5) - 0.5636663150858098975790317e-1 * c(m-4) - 0.4382544850115050165022338e0 * c(m-7) - 0.3806231949664312575822630e0 * c(m-6) 0.5555555555555555555555556e-2 * c(m-9) + 0.1250000000000000000000000e0 * c(m-8) + 0.9698817275172575247533506e0 * c(m-7) + 0.1011699483929608164601067e1 * c(m-5) + 0.1773466968705924819112984e0 * c(m-4) + 0.2239223735771599178951297e-3 * c(m-3) + 0.4325148359756313354830552e0 * c(m-6) -0.2500000000000000000000000e-1 * c(m-8) - 0.3132544850115050165022338e0 * c(m-7) - 0.2322389872063761557916742e0 * c(m-4) - 0.8607044252686413302647675e-3 * c(m-3) - 0.2594851141920572702679681e0 * c(m-6) - 0.6691607091647929161078591e0 * c(m-5) 0.5331362125287625412555844e-1 * c(m-7) + 0.1736073411355428563685818e0 * c(m-5) + 0.8671038084174692625075159e-2 * c(m-3) + 0.8084259844422177692569663e-1 * c(m-6) - 0.1664345989168155800449120e0 * c(m-4) -0.5013247356072127938999311e0 * c(m-5) + 0.5021853752328231128475915e0 * c(m-4) - 0.1197175073672143005877150e-1 * c(m-6) 0.3328335104489738933610597e0 * c(m-5) - 0.3179803804558436283847901e0 * c(m-4) - 0.5111353189352474549563559e-2 * c(m-3) - 0.9741776803777790426705996e-2 * c(m-6) -0.1080774142196007991746827e0 * c(m-5) + 0.9941834083648937298100811e-1 * c(m-4) - 0.4834791406446907590553793e-2 * c(m-3) + 0.1349386478955833378422842e-1 * c(m-6) 0.1319674981073749167009902e-1 * c(m-5) - 0.1060554802883657391328704e-1 * c(m-4) + 0.1911888563316170927411831e-2 * c(m-3) - 0.4503090345217088684223814e-2 * c(m-6);
    %     0.2279091916474916391629437e-1 * c(m-7) + 0.3097763128598982561225538e-1 * c(m-6) - 0.3055555555555555555555556e-1 * c(m-8) - 0.3068985997518740530511593e-1 * c(m-5) - 0.3634246031107139778989373e-2 * c(m-4) 0.5000000000000000000000000e-1 * c(m-8) - 0.3557251496099816106154206e0 * c(m-6) + 0.5490976528821799120017102e-1 * c(m-4) + 0.1838363233410033443348225e0 * c(m-7) + 0.2169790609807602750804271e0 * c(m-5) -0.2500000000000000000000000e-1 * c(m-8) - 0.3132544850115050165022338e0 * c(m-7) - 0.2322389872063761557916742e0 * c(m-4) - 0.8607044252686413302647675e-3 * c(m-3) - 0.2594851141920572702679681e0 * c(m-6) - 0.6691607091647929161078591e0 * c(m-5) 0.5555555555555555555555556e-2 * c(m-8) + 0.1338363233410033443348225e0 * c(m-7) + 0.7391887916719206077121040e0 * c(m-6) + 0.6490333320052011212240632e0 * c(m-4) + 0.3308343404200968256656458e-2 * c(m-3) + 0.2062575706647430620228133e-3 * c(m-2) + 0.3088241944378964404772302e-3 * c(m-1) + 0.4227226173449345042468960e-3 * c(m) + 0.1190362071861893051132274e1 * c(m-5) -0.2720908083525083608370563e-1 * c(m-7) - 0.1931148612480615118957263e0 * c(m-6) - 0.3332941113251635390801278e-1 * c(m-3) + 0.2605582646183255957264249e-3 * c(m-2) + 0.1432116665752147607469646e-2 * c(m-1) + 0.1207544072304193806052558e-2 * c(m) - 0.1348436986667115543203552e1 * c(m-5) + 0.1687723507780044227927853e-1 * c(m-4) 0.3590669644811151307464697e-1 * c(m-6) - 0.5925443480724830632401754e0 * c(m-4) - 0.6897142765790609546343709e-2 * c(m-2) - 0.5363098747528542488971874e-2 * c(m-1) - 0.5205147429855955657625694e-2 * c(m) + 0.9977064356292750529201981e0 * c(m-5) 0.7272438906214475928744770e-1 * c(m-4) + 0.1964682777744275219350831e-1 * c(m-3) - 0.5952475275883259619711594e-2 * c(m-1) - 0.1635430866921887819487473e-2 * c(m) + 0.2921234010758621482958052e-1 * c(m-6) - 0.4659516693228870973898560e0 * c(m-5) 0.5891947149681041048896399e-1 * c(m-4) + 0.1858378996391679448655070e-1 * c(m-3) + 0.7240905383565181316381731e-2 * c(m-2) + 0.2349927974590068869356781e-1 * c(m) - 0.4046360079256766884300687e-1 * c(m-6) + 0.1223513270418807666970488e0 * c(m-5) -0.2404661162020836566908542e-1 * c(m-4) - 0.7348845587775519698437916e-2 * c(m-3) - 0.8105784530576404277872603e-3 * c(m-2) + 0.9574633163221758060736592e-2 * c(m-1) - 0.1828896813877197352675410e-1 * c(m) + 0.1350326632905990039353503e-1 * c(m-6) - 0.1315967038382618382356495e-1 * c(m-5);
    %     0.5522702088127090209264064e-3 * c(m-7) - 0.3265435411305071914756373e-2 * c(m-6) + 0.7462059484530855073291365e-2 * c(m-5) - 0.4748894282038492179461399e-2 * c(m-4) 0.5528052133944605740009217e-1 * c(m-6) - 0.8631683980217122275970376e-1 * c(m-5) - 0.3276463639080639163926118e-1 * c(m-7) + 0.5268984374242044588776166e-1 * c(m-4) 0.5331362125287625412555844e-1 * c(m-7) + 0.1736073411355428563685818e0 * c(m-5) + 0.8671038084174692625075159e-2 * c(m-3) + 0.8084259844422177692569663e-1 * c(m-6) - 0.1664345989168155800449120e0 * c(m-4) -0.2720908083525083608370563e-1 * c(m-7) - 0.1931148612480615118957263e0 * c(m-6) - 0.3332941113251635390801278e-1 * c(m-3) + 0.2605582646183255957264249e-3 * c(m-2) + 0.1432116665752147607469646e-2 * c(m-1) + 0.1207544072304193806052558e-2 * c(m) - 0.1348436986667115543203552e1 * c(m-5) + 0.1687723507780044227927853e-1 * c(m-4) 0.6107825764368264576481962e-2 * c(m-7) + 0.1155752633643216628010304e0 * c(m-6) + 0.2096413329579026439044119e1 * c(m-5) + 0.3357721707576477199985656e0 * c(m-3) + 0.3291545083271862858501887e-3 * c(m-2) + 0.6641183499427826101618457e-2 * c(m-1) + 0.3449455095910233625229891e-2 * c(m) + 0.8270696421223286922584620e0 * c(m-4) -0.4995827370863505253765970e-1 * c(m-6) - 0.1263507837371824205693950e1 * c(m-5) - 0.8712928907711754187084757e-2 * c(m-2) - 0.2487040599390160764166412e-1 * c(m-1) - 0.1486895819265604128572498e-1 * c(m) - 0.1726565392121567634950213e1 * c(m-4) 0.5402985338373433052255418e0 * c(m-5) - 0.1979290298620869974478871e0 * c(m-3) - 0.2760353365637712827793337e-1 * c(m-1) - 0.4671751091575462868310238e-2 * c(m) - 0.6952587985456154591014641e-1 * c(m-6) + 0.1571507277911208446562686e1 * c(m-4) -0.1319358558853174530078498e0 * c(m-5) - 0.1872196143003808021730728e0 * c(m-3) + 0.9147192682075630179962131e-2 * c(m-2) + 0.6712774475803763988977355e-1 * c(m) + 0.9630407686703666967100804e-1 * c(m-6) - 0.6882132817757726722141421e0 * c(m-4) 0.1241625568998496895352046e-1 * c(m-5) + 0.7403484645316174090533193e-1 * c(m-3) - 0.1023976547309387874453988e-2 * c(m-2) + 0.4440063948509876221050939e-1 * c(m-1) - 0.5224403464202056316702078e-1 * c(m) - 0.3213800979246298453953842e-1 * c(m-6) + 0.1178181682424363524005403e0 * c(m-4);
    %     0.6272075574042975468177820e-3 * c(m-6) - 0.1200429618441003833696998e-1 * c(m-5) + 0.1137708862700574079015220e-1 * c(m-4) -0.5373770512016897565958305e-2 * c(m-6) + 0.1453421858063658498587377e0 * c(m-5) - 0.1399684152943489522927794e0 * c(m-4) -0.5013247356072127938999311e0 * c(m-5) + 0.5021853752328231128475915e0 * c(m-4) - 0.1197175073672143005877150e-1 * c(m-6) 0.3590669644811151307464697e-1 * c(m-6) - 0.5925443480724830632401754e0 * c(m-4) - 0.6897142765790609546343709e-2 * c(m-2) - 0.5363098747528542488971874e-2 * c(m-1) - 0.5205147429855955657625694e-2 * c(m) + 0.9977064356292750529201981e0 * c(m-5) -0.4995827370863505253765970e-1 * c(m-6) - 0.1263507837371824205693950e1 * c(m-5) - 0.8712928907711754187084757e-2 * c(m-2) - 0.2487040599390160764166412e-1 * c(m-1) - 0.1486895819265604128572498e-1 * c(m) - 0.1726565392121567634950213e1 * c(m-4) 0.2760393423824887721078848e-1 * c(m-6) + 0.1190550338687608873798462e1 * c(m-5) + 0.4253084328734353394994388e1 * c(m-4) + 0.2306367624634749229113646e0 * c(m-2) + 0.9313657638804699948929701e-1 * c(m-1) + 0.6409299775987186986730499e-1 * c(m) -0.8755807343482262259774782e0 * c(m-5) - 0.3645285178085761821545207e1 * c(m-4) + 0.1033717994630886401730470e0 * c(m-1) + 0.2013769413884797246646959e-1 * c(m) + 0.4106783858513785463625543e-1 * c(m-6) 0.3956598149904136332753521e0 * c(m-5) + 0.1630560443616104907615866e1 * c(m-4) - 0.2421320004064592721552708e0 * c(m-2) - 0.2893557395653431666593814e0 * c(m) - 0.5688529641249387985434413e-1 * c(m-6) -0.7684117160199014594442072e-1 * c(m-5) - 0.2928439026361256842196229e0 * c(m-4) + 0.2710530961648671297733465e-1 * c(m-2) - 0.1662748711097054895317080e0 * c(m-1) + 0.2251991532891353212689574e0 * c(m) + 0.1898341454096471754822498e-1 * c(m-6);
    %     0.9209089963443799485648361e-2 * c(m-5) - 0.3129629392354775191148163e-3 * c(m-6) - 0.8896127024208321966533544e-2 * c(m-4) -0.1059816030196818445908057e0 * c(m-5) + 0.1014880675788250237247178e0 * c(m-4) + 0.4493535440856820866087846e-2 * c(m-6) 0.3328335104489738933610597e0 * c(m-5) - 0.3179803804558436283847901e0 * c(m-4) - 0.5111353189352474549563559e-2 * c(m-3) - 0.9741776803777790426705996e-2 * c(m-6) 0.7272438906214475928744770e-1 * c(m-4) + 0.1964682777744275219350831e-1 * c(m-3) - 0.5952475275883259619711594e-2 * c(m-1) - 0.1635430866921887819487473e-2 * c(m) + 0.2921234010758621482958052e-1 * c(m-6) - 0.4659516693228870973898560e0 * c(m-5) 0.5402985338373433052255418e0 * c(m-5) - 0.1979290298620869974478871e0 * c(m-3) - 0.2760353365637712827793337e-1 * c(m-1) - 0.4671751091575462868310238e-2 * c(m) - 0.6952587985456154591014641e-1 * c(m-6) + 0.1571507277911208446562686e1 * c(m-4) -0.8755807343482262259774782e0 * c(m-5) - 0.3645285178085761821545207e1 * c(m-4) + 0.1033717994630886401730470e0 * c(m-1) + 0.2013769413884797246646959e-1 * c(m) + 0.4106783858513785463625543e-1 * c(m-6) 0.1070920689960817104203947e1 * c(m-5) + 0.3717418466925056542408153e1 * c(m-4) + 0.1166740554279680007487795e0 * c(m-3) + 0.1147318200715868527529827e0 * c(m-1) + 0.6327161147136873807796515e-2 * c(m) + 0.6235373239336055200426697e-1 * c(m-6) -0.6639605248735044787146222e0 * c(m-5) - 0.1865625445986772763641423e1 * c(m-4) + 0.1103611313171476425250639e0 * c(m-3) - 0.9091410269992464604926176e-1 * c(m) - 0.8636954541126674177407762e-1 * c(m-6) 0.1582127073537215443965653e0 * c(m-5) + 0.3746489300753517635549495e0 * c(m-4) - 0.4364163147111892346990101e-1 * c(m-3) - 0.1845476106024151050283847e0 * c(m-1) + 0.7075642937243715046279337e-1 * c(m) + 0.2882271848190011329385407e-1 * c(m-6);
    %     -0.3576545132696983143406173e-2 * c(m-5) + 0.4335019854436220306755673e-3 * c(m-6) + 0.3143043147253361112730605e-2 * c(m-4) 0.4093499466767054661591066e-1 * c(m-5) - 0.3471075437892810033585296e-1 * c(m-4) - 0.6224240288742446280057699e-2 * c(m-6) -0.1080774142196007991746827e0 * c(m-5) + 0.9941834083648937298100811e-1 * c(m-4) - 0.4834791406446907590553793e-2 * c(m-3) + 0.1349386478955833378422842e-1 * c(m-6) 0.5891947149681041048896399e-1 * c(m-4) + 0.1858378996391679448655070e-1 * c(m-3) + 0.7240905383565181316381731e-2 * c(m-2) + 0.2349927974590068869356781e-1 * c(m) - 0.4046360079256766884300687e-1 * c(m-6) + 0.1223513270418807666970488e0 * c(m-5) -0.1319358558853174530078498e0 * c(m-5) - 0.1872196143003808021730728e0 * c(m-3) + 0.9147192682075630179962131e-2 * c(m-2) + 0.6712774475803763988977355e-1 * c(m) + 0.9630407686703666967100804e-1 * c(m-6) - 0.6882132817757726722141421e0 * c(m-4) 0.3956598149904136332753521e0 * c(m-5) + 0.1630560443616104907615866e1 * c(m-4) - 0.2421320004064592721552708e0 * c(m-2) - 0.2893557395653431666593814e0 * c(m) - 0.5688529641249387985434413e-1 * c(m-6) -0.6639605248735044787146222e0 * c(m-5) - 0.1865625445986772763641423e1 * c(m-4) + 0.1103611313171476425250639e0 * c(m-3) - 0.9091410269992464604926176e-1 * c(m) - 0.8636954541126674177407762e-1 * c(m-6) 0.4681819359722749441073885e0 * c(m-5) + 0.1015239189167790053447110e1 * c(m-4) + 0.1043897828092562609502636e0 * c(m-3) + 0.2542001760457345743492403e0 * c(m-2) + 0.1306332157111667628555907e1 * c(m) + 0.1196351539550049336518187e0 * c(m-6) -0.1195777325611201766551392e0 * c(m-5) - 0.2187310061229745694542609e0 * c(m-4) - 0.4128029838349298819825156e-1 * c(m-3) - 0.2845627370491611369031341e-1 * c(m-2) - 0.1016689339350338144430605e1 * c(m) - 0.3992391469197282238624438e-1 * c(m-6);
    %     0.5593983696629863059347067e-3 * c(m-5) - 0.1446656414398166805849327e-3 * c(m-6) - 0.4147327282231696253497740e-3 * c(m-4) -0.6486184157331537899459796e-2 * c(m-5) + 0.4409068609809831485979484e-2 * c(m-4) + 0.2077115547521706413480312e-2 * c(m-6) 0.1319674981073749167009902e-1 * c(m-5) - 0.1060554802883657391328704e-1 * c(m-4) + 0.1911888563316170927411831e-2 * c(m-3) - 0.4503090345217088684223814e-2 * c(m-6) -0.2404661162020836566908542e-1 * c(m-4) - 0.7348845587775519698437916e-2 * c(m-3) - 0.8105784530576404277872603e-3 * c(m-2) + 0.9574633163221758060736592e-2 * c(m-1) - 0.1828896813877197352675410e-1 * c(m) + 0.1350326632905990039353503e-1 * c(m-6) - 0.1315967038382618382356495e-1 * c(m-5) 0.1241625568998496895352046e-1 * c(m-5) + 0.7403484645316174090533193e-1 * c(m-3) - 0.1023976547309387874453988e-2 * c(m-2) + 0.4440063948509876221050939e-1 * c(m-1) - 0.5224403464202056316702078e-1 * c(m) - 0.3213800979246298453953842e-1 * c(m-6) + 0.1178181682424363524005403e0 * c(m-4) -0.7684117160199014594442072e-1 * c(m-5) - 0.2928439026361256842196229e0 * c(m-4) + 0.2710530961648671297733465e-1 * c(m-2) - 0.1662748711097054895317080e0 * c(m-1) + 0.2251991532891353212689574e0 * c(m) + 0.1898341454096471754822498e-1 * c(m-6) 0.1582127073537215443965653e0 * c(m-5) + 0.3746489300753517635549495e0 * c(m-4) - 0.4364163147111892346990101e-1 * c(m-3) - 0.1845476106024151050283847e0 * c(m-1) + 0.7075642937243715046279337e-1 * c(m) + 0.2882271848190011329385407e-1 * c(m-6) -0.1195777325611201766551392e0 * c(m-5) - 0.2187310061229745694542609e0 * c(m-4) - 0.4128029838349298819825156e-1 * c(m-3) - 0.2845627370491611369031341e-1 * c(m-2) - 0.1016689339350338144430605e1 * c(m) - 0.3992391469197282238624438e-1 * c(m-6) 0.3167964748016105299646518e-1 * c(m-5) + 0.4976563420877041544013670e-1 * c(m-4) + 0.1632404042590951953384672e-1 * c(m-3) + 0.3185519088796429015220016e-2 * c(m-2) + 0.2968472090638000742888467e0 * c(m-1) + 0.7912667594695582093926295e0 * c(m) + 0.1332316557164627464149716e-1 * c(m-6);
    % ];
    %
    % M(5,10)=M(10,5);
    % M(m-4,m-9)=M(m-9,m-4);
    %
    % M=M/h;
    %
    % D2=HI*(-M-diag(c)*e_1*S_1+diag(c)*e_m*S_m);

    S2_U = [0.15e2/0.4e1 -0.77e2/0.6e1 0.107e3/0.6e1 -13 0.61e2/0.12e2 -0.5e1/0.6e1;]/h^2;
    S2_1 = zeros(1,m);
    S2_1(1:6) = S2_U;
    S2_m = zeros(1,m);
    S2_m(m-5:m) = fliplr(S2_U);





    % Fourth derivative, 1th order accurate at first 8 boundary points (still
    % yield 5th order convergence if stable: for example u_tt = -u_xxxx

    m4 = 7/240;
    m3 = -2/5;
    m2 = 169/60;
    m1 = -122/15;
    m0 = 91/8;

    M4 = m4*(diag(ones(m-4,1),4)+diag(ones(m-4,1),-4))+m3*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3))+m2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2))+m1*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1))+m0*diag(ones(m,1),0);

    %M4 = (-1/6*(diag(ones(m-3,1),3)+diag(ones(m-3,1),-3) ) + 2*(diag(ones(m-2,1),2)+diag(ones(m-2,1),-2)) -13/2*(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)) + 28/3*diag(ones(m,1),0));

    M4_U = [
        0.1009e4/0.192e3 -0.7657e4/0.480e3 0.9307e4/0.480e3 -0.509e3/0.40e2 0.4621e4/0.960e3 -0.25e2/0.32e2;
        -0.7657e4/0.480e3 0.49513e5/0.960e3 -0.4007e4/0.60e2 0.21799e5/0.480e3 -0.8171e4/0.480e3 0.2657e4/0.960e3;
        0.9307e4/0.480e3 -0.4007e4/0.60e2 0.1399e4/0.15e2 -0.2721e4/0.40e2 0.12703e5/0.480e3 -0.521e3/0.120e3;
        -0.509e3/0.40e2 0.21799e5/0.480e3 -0.2721e4/0.40e2 0.3349e4/0.60e2 -0.389e3/0.15e2 0.559e3/0.96e2;
        0.4621e4/0.960e3 -0.8171e4/0.480e3 0.12703e5/0.480e3 -0.389e3/0.15e2 0.17857e5/0.960e3 -0.1499e4/0.160e3;
        -0.25e2/0.32e2 0.2657e4/0.960e3 -0.521e3/0.120e3 0.559e3/0.96e2 -0.1499e4/0.160e3 0.2225e4/0.192e3;
    ];

    M4(1:6,1:6) = M4_U;

    M4(m-5:m,m-5:m) = flipud( fliplr( M4_U ) );
    M4 = M4/h^3;

    S3_U = [-0.17e2/0.4e1 0.71e2/0.4e1 -0.59e2/0.2e1 0.49e2/0.2e1 -0.41e2/0.4e1 0.7e1/0.4e1;]/h^3;
    S3_1 = zeros(1,m);
    S3_1(1:6) = S3_U;
    S3_m = zeros(1,m);
    S3_m(m-5:m) = fliplr(-S3_U);

    D4 = HI*(M4-e_1*S3_1+e_m*S3_m  + S_1'*S2_1-S_m'*S2_m);


end
