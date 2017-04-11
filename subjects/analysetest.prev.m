%This is the full versino of the original analysis code (in the context of the github version, it crashes at line 567 due to missing stuff).

NBOBJS = 40;
NBNAT=500;
UNITTYPE = 'ndp';
MINCOUNT=4;

%addpath(genpath('/home/tm150/MATLAB/plot2svg_20120520'));
%addpath(genpath('/home/tm150/MATLAB/mi'));

%load('modeloutMatlabHSMax.mat');
%load('modeloutMatlabMax.mat');
%load('modeloutCdiam9size56.mat');
%load('modeloutCdiam9size44MinusIso.mat'); %Not really different

%load('modeloutCdiam9size44_1200_20s_onlyposmod.mat');
%load('modeloutCdiam11size44_1200_20s.mat');
%load('modeloutCdiam13size44_1200_20s.mat');
%load('modeloutCdiam7size44_1200_20s.mat');

load('modelout.mat'); % Standard
%load('../../saccades_test/subjects/modeloutCdiam9size44_1200_20s.mat'); % Standard
%load('modeloutCdiam9size44_1200_20s_STANDARD.mat'); % Standard
%load('modeloutCdiam9size44_1200_20s.mat'); 
%load('modeloutCdiam9size44_1200_20s_Sigma05.mat'); 
%load('modeloutCdiam9size44_1200_20s_DontZero.mat');
%load('modeloutCdiam9size44_1200_20s.mat');

%load('modeloutCdiam7size36_1200_20s.mat');
%load('modeloutCdiam9size44_1200.mat');
%load('modeloutCdiam9size44.mat');
%load('modeloutCdiam9size44Mean.mat');
%load('modeloutCdiam11size48.mat');

%salres = load('modeloutSalDiam9size44.mat');
%rndres = load('modeloutCRanddiam9size44.mat');

%rndres = load('modeloutCdiam9size44_1200_20s_Rand.mat'); % Standard
%salres = load('modeloutCdiam9size44_1200_20s_Sal.mat'); % Standard


EXCLUDE1234 = 1;   % This excludes the first 4 subjects (Andrew, Alissa/Alyssa, Jess and Khadijah).
% The reason for this exclusion is that these subjects saw a buggy version of the code, in which the target-absent trials almost always had objects 1 and 2 (the top-hat and accordion) as a target.
% Note that we still need to include the data in the array for alignment purposes with model output. We just don't use it in the actual analysis.


%load('modeloutCdiam7size44.mat');
%load('modeloutCdiam9size36.mat'); 
%load('modeloutCdiam7size36.mat');
%load('modeloutCdiam13size56.mat');
%load('modeloutC.mat');


results1=[]; results1rep=[]; results2=[]; results2rep=[]; subjs1=[]; subjs1rep=[]; subjs2=[]; subjs2rep=[];

datadir = './subjectsdata/';

%We try to keep the series of 1st blocks similar to that in makepics and attentionMapC, just to check if there is some difference in predictability 
z=load([datadir '/ResultsSaccades_AndrewT_blk1.mat']); 
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*1];
z=load([datadir '/ResultsSaccades_AlissaT_blk2.mat']);  
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*2];
z=load([datadir '/ResultsSaccades_Fetemeh_blk3.mat']);
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*3];
z=load([datadir '/ResultsSaccades_Mario_blk4.mat']);
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*4];
z=load([datadir '/ResultsSaccades_Vivian_blk5.mat']);
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*5];
z=load([datadir '/ResultsSaccades_Sam_blk6.mat']);
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*6];
z=load([datadir '/ResultsSaccades_Claudio_blk7.mat']);
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*7];
z=load([datadir '/ResultsSaccades_Alexandru_blk8.mat']);
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*8];
z=load([datadir '/ResultsSaccades_Mona_blk9.mat']);
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*9];
z=load([datadir '/ResultsSaccades_HC_blk10.mat']);
[x i]=sort(z.sequence); results1=[results1 z.results(i)]; subjs1=[subjs1 ones(1,numel(z.results))*10];


z=load([datadir '/ResultsSaccades_AndrewT2_blk1.mat']); 
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*1];
z=load([datadir '/ResultsSaccades_AlyssaT2_blk2.mat']);  
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*2];
z=load([datadir '/ResultsSaccades_Fetemeh2_blk3.mat']);
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*3];
z=load([datadir '/ResultsSaccades_Mario2_blk4.mat']);
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*4];
z=load([datadir '/ResultsSaccades_Vivian2_blk5.mat']);
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*5];
z=load([datadir '/ResultsSaccades_Sam2_blk6.mat']);
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*6];
z=load([datadir '/ResultsSaccades_Claudio2_blk7.mat']);
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*7];
z=load([datadir '/ResultsSaccades_Alexandru2_blk8.mat']);
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*8];
z=load([datadir '/ResultsSaccades_Mona2_blk9.mat']);
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*9];
z=load([datadir '/ResultsSaccades_HC2_blk10.mat']);
[x i]=sort(z.sequence); results1rep=[results1rep z.results(i)]; subjs1rep=[subjs1rep ones(1,numel(z.results))*10];


z=load([datadir '/ResultsSaccades_Jess_blk1.mat']);  
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*11];
z=load([datadir '/ResultsSaccades_KhadijahH_blk2.mat']);                                                       
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*12];
z=load([datadir '/ResultsSaccades_Christine_blk3.mat']);                                                       
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*13];
z=load([datadir '/ResultsSaccades_GregoryM_blk4.mat']);                                                        
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*14];
z=load([datadir '/ResultsSaccades_Karen_blk5.mat']);                                                           
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*15];
z=load([datadir '/ResultsSaccades_Jenny_blk6.mat']);                                                           
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*16];
z=load([datadir '/ResultsSaccades_Luisa_blk7.mat']);                                                           
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*17];
z=load([datadir '/ResultsSaccades_Nick_blk8.mat']);                                                            
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*18];
z=load([datadir '/ResultsSaccades_Amanda_blk9.mat']);                                                            
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*19];
z=load([datadir '/ResultsSaccades_Wendy_blk10.mat']);                                                            
[x i]=sort(z.sequence); results2=[results2 z.results(i)]; subjs2=[subjs2 ones(1,numel(z.results))*20];


% This one is a bit more complicated because some of the subjects did not come
% back for the second block, or the second block is incomplete (Jess). Yet we
% need data from these blocks to align the array with the model-output data.
% What we do: we include their first block in results2rep, and then we will
% filter these out when we build select2rep

z=load([datadir '/ResultsSaccades_Jess_blk1.mat']);  % TO REMOVE IN SELECT2REP   !!!!!!!!!!
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*11];
z=load([datadir '/ResultsSaccades_KhadijahH_blk2.mat']); % TO REMOVE IN SELECT2REP  !!!!!!!!!                              
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*12];
z=load([datadir '/ResultsSaccades_Christine2_blk3.mat']);                                                                  
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*13];
z=load([datadir '/ResultsSaccades_GregoryM2_blk4.mat']);                                                                   
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*14];
z=load([datadir '/ResultsSaccades_Karen2_blk5.mat']);                                                                      
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*15];
z=load([datadir '/ResultsSaccades_Jenny2_blk6.mat']);                                                                      
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*16];
z=load([datadir '/ResultsSaccades_Luisa2_blk7.mat']);                                                                      
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*17];
z=load([datadir '/ResultsSaccades_Nick_blk8.mat']); % TO REMOVE IN SELECT2REP !!!!!!!!                                     
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*18];
z=load([datadir '/ResultsSaccades_Amanda_blk9.mat']); % TO REMOVE IN SELECT2REP !!!!!!!!                                   
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*19];
z=load([datadir '/ResultsSaccades_Wendy_blk10.mat']); % TO REMOVE IN SELECT2REP !!!!!!!!                                   
[x i]=sort(z.sequence); results2rep=[results2rep z.results(i)]; subjs2rep=[subjs2rep ones(1,numel(z.results))*20];









%results1=[results1(1:370) results1(440:810) results1(880:end)];
%results2=[results2(1:370) results2(440:810) results2(880:end)];
%fixorder=[fixorder(1:370,:); fixorder(440:810,:); fixorder(880:end,:)];
%arraysmodel=[arraysmodel(1:370,:); arraysmodel(440:810,:); arraysmodel(880:end,:)];


% We build the arrays that include relevant values (presented objects, targets, etc.)
% for each series of block presentations. Each block
% was shown at least twice to 2 different subjects (arrays1 and arrays2), and in addition, some
% subjects came back for a second session in which they saw the same block (arrays1rep and arrays2rep). 

arrays1=[]; for i=1:length(results1);  arrays1=[arrays1; results1{i}.probes];  end
isobjinarray1=[]; for n=1:NBOBJS; isobjinarray1 = [isobjinarray1 sum((arrays1==n)')']; end;
tgt1=[]; for i=1:length(results1);  tgt1=[tgt1 results1{i}.target];  end
tgtpos1=[]; for i=1:length(results1);  tgtpos1=[tgtpos1 results1{i}.targetpos];  end
fixedpos1=-1*ones(numel(results1), 6); for i=1:length(results1);  uf=myunique(results1{i}.fixatedobjects, 'stable'); for j=1:numel(uf); fixedpos1(i,j)=uf(j); end;end;
firsts1=[]; for i=1:length(results1); if length(results1{i}.fixatedobjects) > 0; firsts1=[firsts1 results1{i}.fixatedobjects(1)]; else firsts1=[firsts1 -1]; end; end;
fixedobjs1=-1*ones(numel(results1), 6); for i=1:length(results1);  uf=fixedpos1(i,:); for j=1:numel(uf); if uf(j)==-1 continue; end; fixedobjs1(i,j)=arrays1(i, uf(j)); end; end;

arrays1rep=[]; for i=1:length(results1rep);  arrays1rep=[arrays1rep; results1rep{i}.probes];  end
isobjinarray1rep=[]; for n=1:NBOBJS; isobjinarray1rep = [isobjinarray1rep sum((arrays1rep==n)')']; end;
tgt1rep=[]; for i=1:length(results1rep);  tgt1rep=[tgt1rep results1rep{i}.target];  end
tgtpos1rep=[]; for i=1:length(results1rep);  tgtpos1rep=[tgtpos1rep results1rep{i}.targetpos];  end
fixedpos1rep=-1*ones(numel(results1rep), 6); for i=1:length(results1rep);  uf=myunique(results1rep{i}.fixatedobjects, 'stable'); for j=1:numel(uf); fixedpos1rep(i,j)=uf(j); end;end;
firsts1rep=[]; for i=1:length(results1rep); if length(results1rep{i}.fixatedobjects) > 0; firsts1rep=[firsts1rep results1rep{i}.fixatedobjects(1)]; else firsts1rep=[firsts1rep -1]; end; end;
fixedobjs1rep=-1*ones(numel(results1rep), 6); for i=1:length(results1rep);  uf=fixedpos1rep(i,:); for j=1:numel(uf); if uf(j)==-1 continue; end; fixedobjs1rep(i,j)=arrays1rep(i, uf(j)); end; end;


arrays2=[]; for i=1:length(results2);  arrays2=[arrays2; results2{i}.probes];  end
isobjinarray2=[]; for n=1:NBOBJS; isobjinarray2 = [isobjinarray2 sum((arrays2==n)')']; end;
tgt2=[]; for i=1:length(results2);  tgt2=[tgt2 results2{i}.target];  end
tgtpos2=[]; for i=1:length(results2);  tgtpos2=[tgtpos2 results2{i}.targetpos];  end
fixedpos2=-1*ones(numel(results2), 6); for i=1:length(results2);  uf=myunique(results2{i}.fixatedobjects, 'stable'); for j=1:numel(uf); fixedpos2(i,j)=uf(j); end;end;
firsts2=[]; for i=1:length(results2); if length(results2{i}.fixatedobjects) > 0; firsts2=[firsts2 results2{i}.fixatedobjects(1)]; else firsts2=[firsts2 -1]; end; end;
fixedobjs2=-1*ones(numel(results2), 6); for i=1:length(results2);  uf=fixedpos2(i,:); for j=1:numel(uf); if uf(j)==-1 continue; end; fixedobjs2(i,j)=arrays2(i, uf(j)); end; end;
%fixst=[];for i=1:length(results); if length(results{i}.fixatedobjects) > 0; fixst=[fixst results{i}.fixationstarts(1)]; else fixst=[fixst -1]; end; end;


arrays2rep=[]; for i=1:length(results2rep);  arrays2rep=[arrays2rep; results2rep{i}.probes];  end
isobjinarray2rep=[]; for n=1:NBOBJS; isobjinarray2rep = [isobjinarray2rep sum((arrays2rep==n)')']; end;
tgt2rep=[]; for i=1:length(results2rep);  tgt2rep=[tgt2rep results2rep{i}.target];  end
tgtpos2rep=[]; for i=1:length(results2rep);  tgtpos2rep=[tgtpos2rep results2rep{i}.targetpos];  end
fixedpos2rep=-1*ones(numel(results2rep), 6); for i=1:length(results2rep);  uf=myunique(results2rep{i}.fixatedobjects, 'stable'); for j=1:numel(uf); fixedpos2rep(i,j)=uf(j); end;end;
firsts2rep=[]; for i=1:length(results2rep); if length(results2rep{i}.fixatedobjects) > 0; firsts2rep=[firsts2rep results2rep{i}.fixatedobjects(1)]; else firsts2rep=[firsts2rep -1]; end; end;
fixedobjs2rep=-1*ones(numel(results2rep), 6); for i=1:length(results2rep);  uf=fixedpos2rep(i,:); for j=1:numel(uf); if uf(j)==-1 continue; end; fixedobjs2rep(i,j)=arrays2rep(i, uf(j)); end; end;

fixedobjsmod=-1*ones(size(fixorder)); for i=1:size(fixorder,1); uf=fixorder(i,:); for j=1:numel(uf); if uf(j)==-1 continue; end; fixedobjsmod(i,j) = arraysmodel(i, uf(j)); end; end;
fmod = fixedobjsmod(:,1);

%fixedobjssal=-1*ones(size(salres.fixorder)); for i=1:size(salres.fixorder,1); uf=salres.fixorder(i,:); for j=1:numel(uf); if uf(j)==-1 continue; end; fixedobjssal(i,j) = salres.arraysmodel(i, uf(j)); end; end;
%fsal = fixedobjssal(:,1);

%fixedobjsrnd=-1*ones(size(rndres.fixorder)); for i=1:size(rndres.fixorder,1); uf=rndres.fixorder(i,:); for j=1:numel(uf); if uf(j)==-1 continue; end; fixedobjsrnd(i,j) = rndres.arraysmodel(i, uf(j)); end; end;
%frnd = fixedobjsrnd(:,1);

%frnds=[];
%for n=1:1000
%	nbcorrect = sum(fmod(tgtpos1 ~= -10) == tgt1(tgtpos1 ~= -10)'); 
%	fwrong=-1*ones(1,size(fixorder,1)); for i=1:size(fixorder,1); fwrong(i)=tgt1(i); while(fwrong(i)==tgt1(i)) uf=ceil(6*rand());  fwrong(i) = arraysmodel(i, uf); end; end; 
%	selpos= find(tgtpos1 ~= -10); selpos=selpos(randperm(numel(selpos))); selpos=selpos(1:nbcorrect); frnd=fwrong; frnd(selpos)=tgt1(selpos);
%	frnds=[frnds frnd'];
%end;



allowrep2 = subjs2rep ~= 11 & subjs2rep ~= 12 & subjs2rep ~= 18 & subjs2rep ~= 19 & subjs2rep ~= 20; % Exclude those that didn't take the repetition block



% "m" is the main array, collating all the relevant values from both subjects and model.


m=[tgt1' tgtpos1' fixedobjs1(:,1) fixedobjsmod(:,1)  subjs1'; tgt1rep' tgtpos1rep' fixedobjs1rep(:,1) fixedobjsmod(:,1)   subjs1rep'; tgt2' tgtpos2' fixedobjs2(:,1) fixedobjsmod(:,1)  subjs2'; tgt2rep(allowrep2)' tgtpos2rep(allowrep2)' fixedobjs2rep(allowrep2,1) fixedobjsmod(allowrep2,1)  subjs2rep(allowrep2)'];
mallmodobjs=[tgt1' tgtpos1' fixedobjs1(:,1) fixedobjsmod(:,:)  subjs1'; tgt1rep' tgtpos1rep' fixedobjs1rep(:,1) fixedobjsmod(:,:)   subjs1rep'; tgt2' tgtpos2' fixedobjs2(:,1) fixedobjsmod(:,:)  subjs2'; tgt2rep(allowrep2)' tgtpos2rep(allowrep2)' fixedobjs2rep(allowrep2,1) fixedobjsmod(allowrep2,:)  subjs2rep(allowrep2)'];
mbws=[tgt1' tgtpos1' fixedobjs1(:,1) fixedobjs2(:,1) fixedobjsmod(:,1) subjs1' subjs2' ; tgt1rep(allowrep2)' tgtpos1rep(allowrep2)' fixedobjs1rep(allowrep2,1) fixedobjs2rep(allowrep2,1) fixedobjsmod(allowrep2,1) subjs1rep(allowrep2)' subjs2rep(allowrep2)';];
mwis=[tgt1' tgtpos1' fixedobjs1(:,1) fixedobjs1rep(:,1) fixedobjsmod(:,1) subjs1'; tgt2(allowrep2)' tgtpos2(allowrep2)' fixedobjs2(allowrep2,1) fixedobjs2rep(allowrep2,1) fixedobjsmod(allowrep2,1) subjs2(allowrep2)'];
marrays=[tgt1' tgtpos1' fixedobjs1(:,1) fixedobjsmod(:,1)  subjs1' arrays1; tgt1rep' tgtpos1rep' fixedobjs1rep(:,1) fixedobjsmod(:,1)   subjs1rep' arrays1rep; tgt2' tgtpos2' fixedobjs2(:,1) fixedobjsmod(:,1)  subjs2' arrays2; tgt2rep(allowrep2)' tgtpos2rep(allowrep2)' fixedobjs2rep(allowrep2,1) fixedobjsmod(allowrep2,1)  subjs2rep(allowrep2)' arrays2rep(allowrep2,:)];
mwisametgt = mwis; mwidifftgt = mwis;
mbwsametgt = mbws; mbwdifftgt = mbws;
m=m(m(:,3) ~= -1, :);
mallmodobjs=mallmodobjs(mallmodobjs(:,3) ~= -1, :);
marrays=marrays(marrays(:,3) ~= -1, :);
mbws=mbws(mbws(:,3) ~= -1 & mbws(:,4) ~= -1, :); 
mwis=mwis(mwis(:,3) ~= -1 & mwis(:,4) ~= -1, :);
for i=1:15; st=440*(i-1) + 301; mwidifftgt(st:st+69, 4)= mwisametgt(st+70:st+139, 4) ; mwidifftgt(st+70:st+139,4) = mwisametgt(st:st+69, 4); end;
for i=1:15; st=440*(i-1) + 301; mbwdifftgt(st:st+69, 4)= mbwsametgt(st+70:st+139, 4) ; mbwdifftgt(st+70:st+139,4) = mbwsametgt(st:st+69, 4); end;
mwisametgt=mwisametgt(mwisametgt(:,3) ~= -1 & mwisametgt(:,4) ~= -1,:); mwidifftgt=mwidifftgt(mwidifftgt(:,3) ~= -1 & mwidifftgt(:,4) ~= -1,:); 
mbwsametgt=mbwsametgt(mbwsametgt(:,3) ~= -1 & mbwsametgt(:,4) ~= -1,:); mbwdifftgt=mbwdifftgt(mbwdifftgt(:,3) ~= -1 & mbwdifftgt(:,4) ~= -1,:); 

if EXCLUDE1234 == 1
m=m(m(:,5) ~= 1 & m(:,5) ~= 2 & m(:,5) ~= 11 & m(:,5) ~= 12, :);
mallmodobjs=mallmodobjs(mallmodobjs(:,5) ~= 1 & mallmodobjs(:,5) ~= 2 & mallmodobjs(:,5) ~= 11 & mallmodobjs(:,5) ~= 12, :);
mbws=mbws(mbws(:,6) ~= 1 & mbws(:,6) ~= 2 & mbws(:,6) ~= 11 & mbws(:,6) ~= 12 & mbws(:,7) ~= 1 & mbws(:,7) ~= 2 & mbws(:,7) ~= 11 & mbws(:,7) ~= 12, :);
mwis=mwis(mwis(:,6) ~= 1 & mwis(:,6) ~= 2 & mwis(:,6) ~= 11 & mwis(:,6) ~= 12, :);
mwisametgt=mwisametgt(~ismember(mwisametgt(:,6), [1 2 11 12]),:);
mbwsametgt=mbwsametgt(~ismember(mbwsametgt(:,6), [1 2 11 12]) & ~ismember(mbwsametgt(:,7), [1 2 11 12]) ,:);
mwisametgt=mwisametgt(~ismember(mwisametgt(:,6), [1 2 11 12]),:);
mbwdifftgt=mbwdifftgt(~ismember(mbwdifftgt(:,6), [1 2 11 12]) & ~ismember(mbwdifftgt(:,7), [1 2 11 12]) ,:);
marrays=marrays(marrays(:,5) ~= 1 & marrays(:,5) ~= 2 & marrays(:,5) ~= 11 & marrays(:,5) ~= 12, :);
end

ta=m(m(:,2) == -10, :);  tabws=mbws(mbws(:,2) == -10, :);   tawis=mwis(mwis(:,2) == -10, :); taarrays=marrays(marrays(:,2) == -10, :);
tp=m(m(:,2) ~= -10, :);  tpbws=mbws(mbws(:,2) ~= -10, :);   tpwis=mwis(mwis(:,2) ~= -10, :); 
botherr=tp(tp(:,1) ~= tp(:,3) & tp(:,1) ~= tp(:,5), :);
err = tp ~= repmat(tp(:,1), 1, 5);
errbws = tpbws ~= repmat(tpbws(:,1), 1, 7); 
errwis = tpwis ~= repmat(tpwis(:,1), 1, 6);
eqsub = tp == repmat(tp(:,3), 1, 5);

disp('Correlation between errors in both presentations within subjs. (r, p-val):');
[x p] = corr(errwis(:,3), errwis(:,4)); disp(x); disp(p);
disp('Correlation between errors in both presentations between subjs. (r, p-val):');
[x p] = corr(errbws(:,3), errbws(:,4)); disp(x); disp(p);
disp('Correlation between errors in subjects and model (r, p-val):');
[x p] = corr(err(:,3), err(:,4)); disp(x); disp(p);




% Human and model performance
% Success rates after each fixation (excluding target absent trials)
% First we gather the objects actually fixated over all permissible trials, into an Nx6 arrays (N trials, 6 fixations, most of which are -1).
cumperfperobjs=[]; cumperfperobjm=[];
cumperfperobjr=[]; cumperfperobjsal=[];
for tgobj=1:NBOBJS
	so=[fixedobjs1(tgt1'==tgobj & tgtpos1'~=-10 & fixedobjs1(:,1) ~= -1, :); fixedobjs2(tgt2'==tgobj & tgtpos2'~=-10 & fixedobjs2(:,1) ~= -1,:);fixedobjs1rep(tgt1rep'==tgobj & tgtpos1rep'~=-10 & fixedobjs1rep(:,1) ~= -1, :); fixedobjs2rep(allowrep2' & tgt2rep'==tgobj & tgtpos2rep'~=-10 & fixedobjs2rep(:,1) ~= -1,:); ];
	mo=[fixedobjsmod(tgt1'==tgobj & tgtpos1'~=-10 & fixedobjs1(:,1) ~= -1, :); fixedobjsmod(tgt2'==tgobj & tgtpos2'~=-10 & fixedobjs2(:,1) ~= -1,:);fixedobjsmod(tgt1rep'==tgobj & tgtpos1rep'~=-10 & fixedobjs1rep(:,1) ~= -1, :); fixedobjsmod(allowrep2' & tgt2rep'==tgobj & tgtpos2rep'~=-10 & fixedobjs2rep(:,1) ~= -1,:); ];
%	ro=[fixedobjsrnd(tgt1'==tgobj & tgtpos1'~=-10 & fixedobjs1(:,1) ~= -1, :); fixedobjsrnd(tgt2'==tgobj & tgtpos2'~=-10 & fixedobjs2(:,1) ~= -1,:);fixedobjsrnd(tgt1rep'==tgobj & tgtpos1rep'~=-10 & fixedobjs1rep(:,1) ~= -1, :); fixedobjsrnd(allowrep2' & tgt2rep'==tgobj & tgtpos2rep'~=-10 & fixedobjs2rep(:,1) ~= -1,:); ];
%	salo=[fixedobjssal(tgt1'==tgobj & tgtpos1'~=-10 & fixedobjs1(:,1) ~= -1, :); fixedobjssal(tgt2'==tgobj & tgtpos2'~=-10 & fixedobjs2(:,1) ~= -1,:);fixedobjssal(tgt1rep'==tgobj & tgtpos1rep'~=-10 & fixedobjs1rep(:,1) ~= -1, :); fixedobjssal(allowrep2' & tgt2rep'==tgobj & tgtpos2rep'~=-10 & fixedobjs2rep(:,1) ~= -1,:); ];
	tgts=[tgt1(tgt1'==tgobj & tgtpos1'~=-10 & fixedobjs1(:,1) ~= -1) tgt2(tgt2'==tgobj & tgtpos2'~=-10 & fixedobjs2(:,1) ~= -1) tgt1rep(tgt1rep'==tgobj & tgtpos1rep'~=-10 & fixedobjs1rep(:,1) ~= -1) tgt2rep(tgt2rep'==tgobj & allowrep2' & tgtpos2rep'~=-10 & fixedobjs2rep(:,1) ~= -1) ];
	% Now for eqch column of the fixqted-objects array, we count the mean of match between this object and the actual target, and perform the cumulative sum of the resulting means (6, 1 per possible fixation).
	%disp(cumsum(mean(so==repmat(tgts',1,6))));
	%disp(cumsum(mean(mo==repmat(tgts',1,6))));
	cumperfperobjs = [cumperfperobjs; (cumsum(mean(so==repmat(tgts',1,6))))];
	cumperfperobjm = [cumperfperobjm; (cumsum(mean(mo==repmat(tgts',1,6))))];
%	cumperfperobjr = [cumperfperobjr; (cumsum(mean(ro==repmat(tgts',1,6))))];
%	cumperfperobjsal = [cumperfperobjsal; (cumsum(mean(salo==repmat(tgts',1,6))))];
end

figure; 
set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontWeight', 'Bold'); 
errorbar(1:6,mean(cumperfperobjs),std(cumperfperobjs)/sqrt(NBOBJS), 'ko-'); hold on; 
errorbar(1:6,mean(cumperfperobjm),std(cumperfperobjm)/sqrt(NBOBJS), 'k*-'); hold on; 
%errorbar(1:6,mean(cumperfperobjr),std(cumperfperobjr)/sqrt(NBOBJS), 'k*--'); hold on; 
%errorbar(1:6,mean(cumperfperobjsal),std(cumperfperobjsal)/sqrt(NBOBJS), 'ko--'); hold on; 
%plot(cumsum(mean(ro==repmat(tgts',1,6))), 'ko--');plot(cumsum(mean(salo==repmat(tgts',1,6))), 'ko:');
%legend('Subject', 'Model', 'Randomized', 'Saliency model', 'location', 'SouthEast'); title('Cumulative performance after each fixation');
legend('Subject', 'Model'); title('Cumulative performance after each fixation');
axis([.9 6.1 0 1.1]);
%%print('-dpng', '-r300', 'cumperf.png');
%print('-depsc', '-r300', 'cumperf.eps');
% Average number of fixations







%return;





set(0, 'DefaultAxesFontSize', 16); 
wiasametgt=[]; wiadifftgt=[]; bwasametgt=[]; bwadifftgt=[];
mwisametgtTA=mwisametgt(mwisametgt(:,2) == -10, :); mwidifftgtTA=mwidifftgt(mwidifftgt(:,2) == -10, :);
mbwsametgtTA=mbwsametgt(mbwsametgt(:,2) == -10, :); mbwdifftgtTA=mbwdifftgt(mbwdifftgt(:,2) == -10, :);
mbwTP=mbwsametgt(mbwsametgt(:,2) ~= -10, :);  mwiTP=mwisametgt(mwisametgt(:,2) ~= -10, :); 
u=myunique(mwisametgtTA(:,end)); for i=1:length(u); wiasametgtTA(i) = mean(mwisametgtTA(mwisametgtTA(:,end)==u(i), 3) == mwisametgtTA(mwisametgtTA(:,end)==u(i), 4)); end;
u=myunique(mwidifftgtTA(:,end)); for i=1:length(u); wiadifftgtTA(i) = mean(mwidifftgtTA(mwidifftgtTA(:,end)==u(i), 3) == mwidifftgtTA(mwidifftgtTA(:,end)==u(i), 4)); end;
u=myunique(mbwsametgtTA(:,end)); for i=1:length(u); bwasametgtTA(i) = mean(mbwsametgtTA(mbwsametgtTA(:,end)==u(i), 3) == mbwsametgtTA(mbwsametgtTA(:,end)==u(i), 4)); end;
u=myunique(mbwdifftgtTA(:,end)); for i=1:length(u); bwadifftgtTA(i) = mean(mbwdifftgtTA(mbwdifftgtTA(:,end)==u(i), 3) == mbwdifftgtTA(mbwdifftgtTA(:,end)==u(i), 4)); end;
u=myunique(mbwTP(:,end)); for i=1:length(u); bwaTP(i) = mean(mbwTP(mbwTP(:,end)==u(i), 3) == mbwTP(mbwTP(:,end)==u(i), 4)); end;
u=myunique(mwiTP(:,end)); for i=1:length(u); wiaTP(i) = mean(mwiTP(mwiTP(:,end)==u(i), 3) == mwiTP(mwiTP(:,end)==u(i), 4)); end;
%plot2svg('SubjectConsistency.svg');
%Should be a supplementary figure
myf=figure;
subplot(121);errorbar([ mean(wiasametgtTA) mean(wiadifftgtTA)], [std(wiasametgtTA)/sqrt(length(u)) std(wiadifftgtTA)/sqrt(length(u))], '.', 'Linewidth', 2); hold on;
bar([ mean(wiasametgtTA) mean(wiadifftgtTA)]);  set(gca, 'XTick', []); axis([0 3 0 .5]); title('TA, Within Subj');
plot(0:3, .166667*ones(1,4), '--');my_xticklabels(gca,[1 2 3], {{'SAME' 'TGT'} {'DIFF' 'TGT'} {}}, 'Fontsize', 16, 'Fontweight', 'bold'); text(1.2,.45, '***', 'Fontsize', 16, 'Fontweight', 'bold'); line([1 2], [.4 .4], 'linewidth', 2, 'color', 'k');
subplot(122);  errorbar([ mean(bwasametgtTA) mean(bwadifftgtTA)], [std(bwasametgtTA)/sqrt(length(u)) std(bwadifftgtTA)/sqrt(length(u))], '.', 'Linewidth', 2); hold on;
bar([ mean(bwasametgtTA) mean(bwadifftgtTA)]);  set(gca, 'XTick', []); axis([0 3 0 .5]); title('TA, Between Subj' );
plot(0:3, .166667*ones(1,4), '--');  my_xticklabels(gca,[1 2], {{'SAME' 'TGT'} {'DIFF' 'TGT'} }, 'FontSize', 16, 'Fontweight', 'bold'); text(1.3,.45, '***', 'Fontsize', 16, 'Fontweight', 'bold'); line([1 2], [.4 .4], 'linewidth', 2, 'color', 'k');
disp('P-val of ranksum for bwasametgtTA and bwadifftgtTA (8 and 8 vals)'); disp(ranksum(bwasametgtTA,bwadifftgtTA));
disp('P-val of ranksum for wiasametgtTA and wiadifftgtTA (13 and 13 vals)'); disp(ranksum(wiasametgtTA,wiadifftgtTA));
%pos=get(myf, 'Position'); set(myf, 'Position', [pos(1) pos(2) pos(3)*.5 pos(4)]);
%print('-dpng', '-r300', 'consistency.png');
%print('-depsc', '-r300', 'consistency.eps');



set(0, 'DefaultAxesFontSize', 16); 
perfmodtpagg=mean(tp(:,1)==tp(:,4));perfsubtpagg=mean(tp(:,1)==tp(:,3));
perfmodtaagg=mean(ta(:,1)==ta(:,4));perfsubtaagg=mean(ta(:,1)==ta(:,3));
for nums=1:20
	msub{nums} = m(m(:,end)==nums, :);
	mbwssub{nums} = mbws(mbws(:,end)==nums | mbws(:, end-1)==nums, :);
	mwissub{nums} = mwis(mwis(:,end)==nums, :);
	tpsub{nums} = tp(tp(:,end)==nums, :);
	tpbwssub{nums} = tpbws(tpbws(:,end)==nums | tpbws(:, end-1)==nums, :);
	tpwissub{nums} = tpwis(tpwis(:,end)==nums, :);
	tasub{nums} = ta(ta(:,end)==nums, :);
	tabwssub{nums} = tabws(tabws(:,end)==nums | tabws(:, end-1)==nums, :);
	tawissub{nums} = tawis(tawis(:,end)==nums, :);
	wiaallsub{nums} = mean(mwissub{nums}(:,4) == mwissub{nums}(:,3));	
	bwaallsub{nums} = mean(mbwssub{nums}(:,4) == mbwssub{nums}(:,3));
	msaallsub{nums} = mean(msub{nums}(:,4) == msub{nums}(:,3)); 
	wiatasub{nums} = mean(tawissub{nums}(:,4) == tawissub{nums}(:,3));	wiatpsub{nums} = mean(tpwissub{nums}(:,4) == tpwissub{nums}(:,3)); wiasub{nums} = mean(mwissub{nums}(:,3) == mwissub{nums}(:,4));
	bwatasub{nums} = mean(tabwssub{nums}(:,4) == tabwssub{nums}(:,3));  bwatpsub{nums} = mean(tpbwssub{nums}(:,4) == tpbwssub{nums}(:,3)); bwasub{nums} = mean(mbwssub{nums}(:,4) == mbwssub{nums}(:,3));
	msatasub{nums} = mean(tasub{nums}(:,4) == tasub{nums}(:,3));  msatpsub{nums} = mean(tpsub{nums}(:,4) == tpsub{nums}(:,3)); msasub{nums} = mean(msub{nums}(:,4) == msub{nums}(:,3));
	bwabesub{nums} = mean(tpbwssub{nums}(tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,3) & tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,4), 3) == tpbwssub{nums}(tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,3) & tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,4), 4));
	wiabesub{nums} = mean(tpwissub{nums}(tpwissub{nums}(:,1) ~= tpwissub{nums}(:,3) & tpwissub{nums}(:,1) ~= tpwissub{nums}(:,4), 3) == tpwissub{nums}(tpwissub{nums}(:,1) ~= tpwissub{nums}(:,3) & tpwissub{nums}(:,1) ~= tpwissub{nums}(:,4), 4));
	msabesub{nums} = mean(tpsub{nums}(tpsub{nums}(:,1) ~= tpsub{nums}(:,3) & tpsub{nums}(:,1) ~= tpsub{nums}(:,4), 3) == tpsub{nums}(tpsub{nums}(:,1) ~= tpsub{nums}(:,3) & tpsub{nums}(:,1) ~= tpsub{nums}(:,4), 4));
	perfsub{nums} = mean(msub{nums}(:,1) == msub{nums}(:,3));
	perfmodsub{nums} = mean(msub{nums}(:,1) == msub{nums}(:,4));
	perfsubrep1{nums} = mean(mwissub{nums}(:,1) == mwissub{nums}(:,3));	perfsubrep2{nums} = mean(mwissub{nums}(:,1) == mwissub{nums}(:,4)); perfsub1sub{nums} = mean(mbwssub{nums}(:,1) == mbwssub{nums}(:,3)); perfsub2sub{nums} = mean(mbwssub{nums}(:,1) == mbwssub{nums}(:,4));
	% 'chance' agreement due to performance alone (on TP trials) = pA(right)*pB(right) + 5 * (1-pA(right))/5 * (1-pB(right))/5
	chancemsaallsub{nums} = perfsub{nums} * perfmodsub{nums} + (1-perfsub{nums})*(1-perfmodsub{nums})/5;
	chancewiaallsub{nums} = perfsubrep1{nums} * perfsubrep2{nums} + (1-perfsubrep1{nums})*(1-perfsubrep2{nums})/5;
	chancebwaallsub{nums} = perfsub1sub{nums} * perfsub2sub{nums} + (1-perfsub1sub{nums})*(1-perfsub2sub{nums})/5;
	% P-values are determined by the cumulative binomial distribution: what is the probability of obtaining that many matches minus one, or fewer, within this number of trials, if the match probability is determined by chance? The p-value is 1 minus that.  
	pvalwiaallsub{nums}=1-binocdf(-1+sum(mwissub{nums}(:,4) == mwissub{nums}(:,3)), numel(mwissub{nums}(:,4)), chancewiaallsub{nums});
	pvalwiatasub{nums}=1-binocdf(-1+sum(tawissub{nums}(:,4) == tawissub{nums}(:,3)), numel(tawissub{nums}(:,4)), 1/6);
	pvalwiabesub{nums}=1-binocdf(-1+sum(tpwissub{nums}(tpwissub{nums}(:,1) ~= tpwissub{nums}(:,3) & tpwissub{nums}(:,1) ~= tpwissub{nums}(:,4),4) == tpwissub{nums}(tpwissub{nums}(:,1) ~= tpwissub{nums}(:,3) & tpwissub{nums}(:,1) ~= tpwissub{nums}(:,4),3)), numel(tpwissub{nums}(tpwissub{nums}(:,1) ~= tpwissub{nums}(:,3) & tpwissub{nums}(:,1) ~= tpwissub{nums}(:,4),4)), .2);
	pvalbwaallsub{nums}=1-binocdf(-1+sum(mbwssub{nums}(:,4) == mbwssub{nums}(:,3)), numel(mbwssub{nums}(:,4)), chancebwaallsub{nums});
	pvalbwatasub{nums}=1-binocdf(-1+sum(tabwssub{nums}(:,4) == tabwssub{nums}(:,3)), numel(tabwssub{nums}(:,4)), 1/6);
	pvalbwabesub{nums}=1-binocdf(-1+sum(tpbwssub{nums}(tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,3) & tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,4),4) == tpbwssub{nums}(tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,3) & tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,4),3)), numel(tpbwssub{nums}(tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,3) & tpbwssub{nums}(:,1) ~= tpbwssub{nums}(:,4),4)), .2);
	pvalmsaallsub{nums}=1-binocdf(-1+sum(msub{nums}(:,4) == msub{nums}(:,3)), numel(msub{nums}(:,4)), chancemsaallsub{nums});
	pvalmsatasub{nums}=1-binocdf(-1+sum(tasub{nums}(:,4) == tasub{nums}(:,3)), numel(tasub{nums}(:,4)), 1/6);
	pvalmsabesub{nums}=1-binocdf(-1+sum(tpsub{nums}(tpsub{nums}(:,1) ~= tpsub{nums}(:,3) & tpsub{nums}(:,1) ~= tpsub{nums}(:,4),4) == tpsub{nums}(tpsub{nums}(:,1) ~= tpsub{nums}(:,3) & tpsub{nums}(:,1) ~= tpsub{nums}(:,4),3)), numel(tpsub{nums}(tpsub{nums}(:,1) ~= tpsub{nums}(:,3) & tpsub{nums}(:,1) ~= tpsub{nums}(:,4),4)), .2);
end
chancemsatpagg= perfmodtpagg * perfsubtpagg + (1-perfmodtpagg) * (1-perfsubtpagg)/5; 
chancemsataagg= perfmodtaagg * perfsubtaagg + (1-perfmodtaagg) * (1-perfsubtaagg)/5; 
myf=figure;
bwatasub=bwatasub(1:end/2); bwaallsub=bwaallsub(1:end/2); bwabesub=bwabesub(1:end/2); 
pvalbwatasub=pvalbwatasub(1:end/2);pvalbwaallsub=pvalbwaallsub(1:end/2);pvalbwabesub=pvalbwabesub(1:end/2);
chancebwaallsub=chancebwaallsub(1:end/2);
subplot(431); bar(cell2mat(wiaallsub)); hold on; chc=plot(cell2mat(chancewiaallsub),'r+', 'Linewidth', 1); ylabel('WITHIN-SUBJECT '); chcperf=plot(0:21,ones(1,22)./6,'r', 'Linewidth', 1); plot(find(~isnan(cell2mat(wiaallsub)) & cell2mat(pvalwiaallsub)<.05), .75, 'k*'); title('ALL TRIALS'); axis([0 21 0 1]);
subplot(432); bar(cell2mat(wiatasub)); hold on; plot(0:21,ones(1,22)./6,'r', 'Linewidth', 1);  plot(find(~isnan(cell2mat(wiatasub)) & cell2mat(pvalwiatasub)<.05), .6, 'k*'); axis([0 21 0 1]); title('TARGET-ABSENT TRIALS'); 
subplot(433); bar(cell2mat(wiabesub)); hold on; plot(0:21,.2*ones(1,22),'r', 'Linewidth', 1); plot(find(cell2mat(pvalwiabesub)< .05 & ~isnan(cell2mat(wiabesub))), .9, 'k*');title('BOTH-WRONG TRIALS'); axis([0 21 0 1]);
subplot(434); bar(cell2mat(bwaallsub)); hold on; plot(0:11,ones(1,12)./6,'r', 'Linewidth', 1);  plot(cell2mat(chancebwaallsub),'r+', 'Linewidth', 1); plot(find(~isnan(cell2mat(bwaallsub)) & cell2mat(pvalbwaallsub) < .05) , .75, 'k*'); ylabel('BETWEEN-SUBJECTS '); axis([0 11 0 1]);
subplot(435); bar(cell2mat(bwatasub)); hold on; plot(0:11,ones(1,12)./6,'r', 'Linewidth', 1);  plot(find(~isnan(cell2mat(bwatasub)) & cell2mat(pvalbwatasub) < .05) , .6, 'k*'); axis([0 11 0 1]);
subplot(436); bar(cell2mat(bwabesub)); hold on; plot(0:11,.2*ones(1,12),'r', 'Linewidth', 1); plot(find(~isnan(cell2mat(bwabesub)) & cell2mat(pvalbwabesub)  < .05), .6, 'k*');axis([0 11 0 1]);
subplot(437); bar(cell2mat(msaallsub)); hold on;chc=plot(0:21,ones(1,22)./6,'r', 'Linewidth', 1);  chcperf=plot(cell2mat(chancemsaallsub),'r+', 'Linewidth', 1); siglev=plot(find(~isnan(cell2mat(msaallsub)) & cell2mat(pvalmsaallsub)< .05), .75, 'k*'); ylabel('MODEL-SUBJECT '); axis([0 21 0 1]);
subplot(438); bar(cell2mat(msatasub)); hold on; plot(0:21,ones(1,22)./6,'r', 'Linewidth', 1); plot(find(~isnan(cell2mat(msatasub)) & cell2mat(pvalmsatasub)< .05), .75, 'k*');axis([0 21 0 1]);
legend([chc, chcperf, siglev(1)], 'Chance', 'Chance level based on performance alone', 'p<.05', 'Location', 'SouthOutside'); 
subplot(439); bar(cell2mat(msabesub)); hold on; plot(0:21,.2*ones(1,22),'r', 'Linewidth', 1); siglev = plot(find(~isnan(cell2mat(msabesub)) & cell2mat(pvalmsabesub) < .05), .6, 'k*');axis([0 21 0 1]);
pos=get(myf, 'Position'); set(myf, 'Position', [pos(1) pos(2) 2*pos(3) 3*pos(4)]);

%figure; 
%%set(0, 'DefaultAxesFontSize', 11); 
%title({'Mean and s.e.m of difference from chance across subjects'});
%subplot(131); hold on;
%z1=cell2mat(wiatpsub)-cell2mat(chancewiatpsub); z1=z1(~isnan(z1)); %signtest(z1)
%z2=cell2mat(bwatpsub)-cell2mat(chancebwatpsub); 
%z3=cell2mat(msatpsub)-cell2mat(chancemsatpsub); 
%errorbar([mean(z1) mean(z2) mean(z3)] , [std(z1)/sqrt(numel(z1)) std(z2)/sqrt(numel(z2)) std(z3)/sqrt(numel(z3))], '.', 'linewidth', 2); bar([mean(z1) mean(z2) mean(z3)]); 
%axis([0 4 0 .3]);
%%set(gca, 'XTick', [1 2 3]);set(gca, 'XTickLabel', {'TP' 'TA' 'BE'});
%ylabel('Target Present');
%subplot(132); hold on;
%z1=cell2mat(wiatasub)- .1666667; z1=z1(~isnan(z1));
%z2=cell2mat(bwatasub)- .1666667;
%z3=cell2mat(msatasub)- .1666667;
%errorbar([mean(z1) mean(z2) mean(z3)] , [std(z1)/sqrt(numel(z1)) std(z2)/sqrt(numel(z2)) std(z3)/sqrt(numel(z3))], '.', 'linewidth', 2); bar([mean(z1) mean(z2) mean(z3)]); 
%axis([0 4 0 .3]);
%%set(gca, 'XTick', [1 2 3]);set(gca, 'XTickLabel', {'TP' 'TA' 'BE'});
%ylabel('Target Absent');
%subplot(133); hold on;
%z1=cell2mat(wiabesub)- .2; z1=z1(~isnan(z1));
%z2=cell2mat(bwabesub)- .2;
%z3=cell2mat(msabesub)- .2;
%errorbar([mean(z1) mean(z2) mean(z3)] , [std(z1)/sqrt(numel(z1)) std(z2)/sqrt(numel(z2)) std(z3)/sqrt(numel(z3))], '.', 'linewidth', 2); bar([mean(z1) mean(z2) mean(z3)]); 
%%set(gca, 'XTick', [1 2 3]);set(gca, 'XTickLabel', {'TP' 'TA' 'BE'});
%text(1, .2, '***', 'FontSize',14);
%ylabel('Both-Error');
%axis([0 4 0 .3]);


myf=figure;
set(0, 'DefaultAxesFontSize', 12, 'DefaultAxesFontWeight', 'Bold'); 
subplot(141);
msa=cell2mat(msaallsub); bwa=cell2mat(bwaallsub); wia=cell2mat(wiaallsub); wia=wia(wia>0);
wia=wia(wia>0); msa=msa(msa>0); bwa=bwa(bwa>0);
errorbar([.5 2 3 ], [mean(wia) mean(bwa) mean(msa) ],...
				[std(wia)./sqrt(numel(wia)), std(bwa)./sqrt(numel(bwa)) std(msa)./sqrt(numel(msa))], ...
				'.', 'linewidth', 1); hold on;
b1=bar(.5, mean(wia), 'BarWidth', .3, 'FaceColor',[.2 .5 1], 'Edgecolor', [.5 .5 .5]); hold on;
b2=bar(2, mean(bwa), 'BarWidth', .8, 'FaceColor',[0 0 1]); 
b3=bar(3, mean(msa),  'FaceColor',[1 0 0]); 
b4=plot(0:4, ones(1,5)./6, '--');
set(gca, 'XTick', []);
axis([0 4 0 .6]);
title('All Trials');
subplot(142)
msa=cell2mat(msatasub); bwa=cell2mat(bwatasub); wia=cell2mat(wiatasub); wia=wia(wia>0);
wia=wia(wia>0); msa=msa(msa>0); bwa=bwa(bwa>0);
errorbar([.5 2 3 ], [mean(wia) mean(bwa) mean(msa) ],...
				[std(wia)./sqrt(numel(wia)), std(bwa)./sqrt(numel(bwa)) std(msa)./sqrt(numel(msa))], ...
				'.', 'linewidth', 1); hold on;
bar(.5, mean(wia), 'BarWidth', .3, 'FaceColor',[.2 .5 1], 'Edgecolor', [.5 .5 .5]); hold on;
bar(2, mean(bwa), 'BarWidth', .8, 'FaceColor',[0 0 1]); 
bar(3, mean(msa),  'FaceColor',[1 0 0]); 
plot(0:4, ones(1,5)./6, '--');
set(gca, 'XTick', []);
%axis([0 4 0 .5]);
axis([0 4 0 .6]);
title('Target Absent');

bwaTA = bwa; msaTA = msa; wiaTA = wia;
disp(['Target-absent agreements: Between subjects (mean, std): ' num2str(mean(bwaTA)) ' +/- ' num2str(std(bwaTA))  ', model-subjects (mean, std): ' num2str(mean(msaTA)) ' +/-' num2str(std(msaTA))]);

subplot(143)
msa=cell2mat(msabesub); bwa=cell2mat(bwabesub); wia=cell2mat(wiabesub); wia=wia(wia>0); msa=msa(msa>0); bwa=bwa(bwa>0);
errorbar([.5 2 3 ], [mean(wia) mean(bwa) mean(msa) ],...
				[std(wia)./sqrt(numel(wia)), std(bwa)./sqrt(numel(bwa)) std(msa)./sqrt(numel(msa))], ...
				'.', 'linewidth', 1); hold on;
b1=bar(.5, mean(wia), 'BarWidth', .3, 'FaceColor',[.2 .5 1], 'Edgecolor', [.5 .5 .5]); hold on;
b2=bar(2, mean(bwa), 'BarWidth', .8, 'FaceColor',[0 0 1]); 
b3=bar(3, mean(msa),  'FaceColor',[1 0 0]); 
b4=plot(0:4, ones(1,5)./5, '--');  % /5 because both-error means one option is removed (the target), leaving only 5 possible choices
set(gca, 'XTick', []);
axis([0 4 0 .6]);
title('Both-Error');
legend([b1 b2 b3 b4], 'Within-subject', 'Between-subject', 'Model-subject', 'Chance', 'Location', 'EastOutside'); legend boxoff;

set(myf, 'Position', [100 100 1000 200]);

 set(gcf,'PaperPositionMode','auto');

%print('-depsc', '-r300', 'modsub.eps');
%%print('-dpng', '-r300', 'modsub.png');
%pos=get(myf, 'Position'); set(myf, 'Position', [pos(1) pos(2) pos(3)*1.1 pos(4)]);
%pos=get(myf, 'Position'); set(myf, 'Position', [pos(1:2) ceil(2*pos(3)/3) pos(4)]);


marrays=[tgt1' tgtpos1' fixedobjs1(:,1) fixedobjsmod(:,1)  subjs1' arrays1; tgt1rep' tgtpos1rep' fixedobjs1rep(:,1) fixedobjsmod(:,1)   subjs1rep' arrays1rep; tgt2' tgtpos2' fixedobjs2(:,1) fixedobjsmod(:,1)  subjs2' arrays2; ...
tgt2rep(allowrep2)' tgtpos2rep(allowrep2)' fixedobjs2rep(allowrep2,1) fixedobjsmod(allowrep2,1)  subjs2rep(allowrep2)' arrays2rep(allowrep2,:)];
marrays=marrays(marrays(:,3) ~= -1, :);
if EXCLUDE1234 == 1
marrays=marrays(marrays(:,5) ~= 1 & marrays(:,5) ~= 2 & marrays(:,5) ~= 11 & marrays(:,5) ~= 12, :);
end
tp=marrays(marrays(:,2) ~= -10, :);
taarrays=marrays(marrays(:,2) == -10, :);
botherr=tp(tp(:,1) ~= tp(:,3) & tp(:,1) ~= tp(:,4), :);


% Confusion matrices for positions
% We compute the position of the subject's first saccade, and of the model's
% first saccade, for all trials.
%
% But BE CAREFUL... since object position is randomized... this can only work for model-subject comparison, NOT for between-humans comparisons!!!!
%
z= taarrays(:,end-5:end)  == repmat(taarrays(:,4),1,6); [modpos j] = ind2sub( size(z'), find(z'));
z= taarrays(:,end-5:end)  == repmat(taarrays(:,3),1,6); [subpos j] = ind2sub( size(z'), find(z'));
matposta=zeros(6); 
for i=1:size(taarrays, 1); matposta(subpos(i), modpos(i)) =  matposta(subpos(i), modpos(i)) +1; end;
matposta=matposta./ repmat(sum(matposta')', 1, 6);  
matposold=matposta;

z= marrays(:,end-5:end)  == repmat(marrays(:,4),1,6); [modpos j] = ind2sub( size(z'), find(z'));
z= marrays(:,end-5:end)  == repmat(marrays(:,3),1,6); [subpos j] = ind2sub( size(z'), find(z'));
matposall=zeros(6); 
for i=1:size(marrays, 1); matposall(subpos(i), modpos(i)) =  matposall(subpos(i), modpos(i)) +1; end;
matposall=matposall./ repmat(sum(matposall')', 1, 6);  

z= botherr(:,end-5:end)  == repmat(botherr(:,4),1,6); [modpos j] = ind2sub( size(z'), find(z'));
z= botherr(:,end-5:end)  == repmat(botherr(:,3),1,6); [subpos j] = ind2sub( size(z'), find(z'));
matposbe=zeros(6); 
for i=1:size(botherr, 1); matposbe(subpos(i), modpos(i)) =  matposbe(subpos(i), modpos(i)) +1; end;
matposbe=matposbe./ repmat(sum(matposbe')', 1, 6);  



%return;


myf=figure;
subplot(131);
imagesc(matposall); set(gca, 'XTick', []);set(gca, 'YTick', []);xlabel('Model'); ylabel('Subject');title('All Trials');
cmapgreen=colormap(gray); cmapgreen(:,[1 3])= cmapgreen(:,[1 3])/2; colormap(cmapgreen); CB=colorbar;
subplot(132);
imagesc(matposta); set(gca, 'XTick', []);set(gca, 'YTick', []);xlabel('Model'); ylabel('Subject');title('Target Absent');
cmapgreen=colormap(gray); cmapgreen(:,[1 3])= cmapgreen(:,[1 3])/2; colormap(cmapgreen); CB=colorbar;
subplot(133);
imagesc(matposbe); set(gca, 'XTick', []);set(gca, 'YTick', []);xlabel('Model'); ylabel('Subject');title('Both-error');
cmapgreen=colormap(gray); cmapgreen(:,[1 3])= cmapgreen(:,[1 3])/2; colormap(cmapgreen); CB=colorbar;

set(myf, 'Position', [100 100 1500 300]);
 set(gcf,'PaperPositionMode','auto');

%print('-depsc', '-r300', 'matpos.eps');





%return;




imgdata={}; 
for i=1:NBOBJS
		%disp(['Isolated image for object ' num2str(i)]);
	    imgdata{i}=imread(['images/obj' num2str(i) '.png']);
		diso{i} = load(['./out/C2bout_ndp_obj' num2str(i) '.png.raw.out'])';
end;
filts = buildv1filts();
 imgfilt={}; orenergies={};histobj={};
for i=1:NBOBJS
	%disp(['Isolated image for object ' num2str(i)]); orenergies{i}=zeros(1,4);
	for ori=1:4
		imgfilt{i}{ori} = imfilter(imgdata{i}, filts{ori}, 'replicate');
		z = imgfilt{i}{ori}(imgdata{i}(:)<127 | imgdata{i}(:)>129);
		orenergies{i}(ori) = mean(z);
	end
	histobj{i} = hist(imgdata{i}(:),1:255); histobj{i}(127:129)=0;
end;
objsizes=zeros(1,NBOBJS); objlums=zeros(1,NBOBJS);
for i=1:NBOBJS
	objsizes(i)=sum(imgdata{i}(:)<127 | imgdata{i}(:)>129);
	objlums(i)=mean(imgdata{i}(imgdata{i}<127 | imgdata{i}>129));
end;
objsizeabsdevs=abs(objsizes - mean(objsizes)); objlumabsdevs=abs(objlums - mean(objlums));
matcorr=zeros(NBOBJS);matcorrhist=zeros(NBOBJS); matcorror=zeros(NBOBJS);matobjsizes=zeros(NBOBJS);matobjlums=zeros(NBOBJS);matabsdiffsize=zeros(NBOBJS);matobjcontrasts=zeros(NBOBJS); matobjmeanabs=zeros(NBOBJS); matobjmsdiffbg=zeros(NBOBJS);matdist=zeros(NBOBJS);
matobjsizeabsdevs=zeros(NBOBJS); matobjlumabsdevs=zeros(NBOBJS); matdistC2=zeros(NBOBJS);
for so=1:NBOBJS
	matobjsizes(so,:)=objsizes(so);  matobjmeanabs(so,:)=mean(abs(imgdata{so}(:)-128));
	matobjlums(so,:)=objlums(so);
	matobjsizeabsdevs(so,:)=objsizeabsdevs(so); matobjlumabsdevs(so,:)=objlumabsdevs(so);
	tmp= imgdata{so}(:); tmp=tmp(tmp<127 | tmp>129); tmp=tmp-mean(tmp); matobjcontrasts(so,:)=mean(tmp.^2);
	tmp= imgdata{so}(:); tmp=tmp(tmp<127 | tmp>129); tmp=tmp-128; matobjmsdiffbg(so,:)=mean(tmp.^2);
end
for so=1:NBOBJS
	for to=1:NBOBJS
		matdistC2(so, to) = sqrt(sum((diso{so}-diso{to}).^2));
		matdist(so, to) = sqrt(sum((double(imgdata{so}(:))-double(imgdata{to}(:))).^2));
		matcorr(so, to) = corr(double(imgdata{so}(:)),double(imgdata{to}(:)));
		matcorrhist(so, to) = corr(double(histobj{so}(:)),double(histobj{to}(:)));
		matcorror(so, to) = corr(orenergies{so}(:), orenergies{to}(:));
		matabsdiffsize(so, to) = abs(objsizes(so)-objsizes(to));
		matabsdifflum(so, to) = abs(objlums(so) - objlums(to));
	end;
end;


marrays=[tgt1' tgtpos1' fixedobjs1(:,1) fixedobjsmod(:,1)  subjs1' arrays1; tgt1rep' tgtpos1rep' fixedobjs1rep(:,1) fixedobjsmod(:,1)   subjs1rep' arrays1rep; tgt2' tgtpos2' fixedobjs2(:,1) fixedobjsmod(:,1)  subjs2' arrays2; tgt2rep(allowrep2)' tgtpos2rep(allowrep2)' fixedobjs2rep(allowrep2,1) fixedobjsmod(allowrep2,1)  subjs2rep(allowrep2)' arrays2rep(allowrep2,:)];
marrays1=[tgt1' tgtpos1' fixedobjs1(:,1) fixedobjsmod(:,1)  subjs1' arrays1; tgt1rep(allowrep2)' tgtpos1rep(allowrep2)' fixedobjs1rep(allowrep2,1) fixedobjsmod(allowrep2,1)   subjs1rep(allowrep2)' arrays1rep(allowrep2,:)];
marrays2=[tgt2' tgtpos2' fixedobjs2(:,1) fixedobjsmod(:,1)  subjs2' arrays2; tgt2rep(allowrep2)' tgtpos2rep(allowrep2)' fixedobjs2rep(allowrep2,1) fixedobjsmod(allowrep2,1)   subjs2rep(allowrep2)' arrays2rep(allowrep2,:)];
mattvals=[usattvals; usattvals; usattvals; usattvals(allowrep2,:)];

sel = marrays(:,2)==-10 & marrays(:,3)~=-1;
sel12 = marrays1(:,2) == -10 & marrays1(:,3) ~= -1 & marrays2(:,2) == -10 & marrays2(:,3) ~= -1;

%EXCLUDE1234=0;
if EXCLUDE1234 == 1
	sel = sel & (marrays(:,5) ~= 1 & marrays(:,5) ~= 2 & marrays(:,5) ~= 11 & marrays(:,5) ~= 12);
	sel12 = sel12 & marrays1(:,5) ~= 1 & marrays1(:,5) ~= 2 & marrays1(:,5) ~= 11 & marrays1(:,5) ~= 12 ...
				 & marrays2(:,5) ~= 1 & marrays2(:,5) ~= 2 & marrays2(:,5) ~= 11 & marrays2(:,5) ~= 12;
end

taattvals=mattvals(sel,:);
taarrays=marrays(sel, :);
taarrays1=marrays1(sel12,:);
taarrays2=marrays2(sel12,:);

subsel=taarrays(:, 6:end)==repmat(taarrays(:,3), 1, 6); subobjs=taarrays(:,3); [i j]=ind2sub(size(subsel'), find(subsel')); subpos=i;
subsel1=taarrays1(:, 6:end)==repmat(taarrays1(:,3), 1, 6); subobjs1=taarrays1(:,3); [i j]=ind2sub(size(subsel1'), find(subsel1'));subpos1=i;
modsel1=taarrays1(:, 6:end)==repmat(taarrays1(:,4), 1, 6); modobjs1=taarrays1(:,4); [i j]=ind2sub(size(modsel1'), find(modsel1'));modpos1=i;
subsel2=taarrays1(:, 6:end)==repmat(taarrays2(:,3), 1, 6); subobjs2=taarrays2(:,3); [i j]=ind2sub(size(subsel2'), find(subsel2'));subpos2=i;% YES, subsel2=taarrays1 - because we want to compare with subsel1, and the array positions are different
modsel=taarrays(:, 6:end)==repmat(taarrays(:,4), 1, 6); modobjs=taarrays(:,4); [i j]=ind2sub(size(modsel'), find(modsel'));modpos=i;
histcorrvals=[]; orcorrvals=[]; absdiffsizevals=[]; absdifflumvals=[]; corrvals=[];
histcorrvals1=[]; orcorrvals1=[]; absdiffsizevals1=[]; absdifflumvals1=[]; corrvals1=[];
for i=1:6
	corrvals=[corrvals matcorr(sub2ind(size(matcorr), taarrays(:,1), taarrays(:,5+i)))];
	histcorrvals=[histcorrvals matcorrhist(sub2ind(size(matcorrhist), taarrays(:,1), taarrays(:,5+i)))];
	absdiffsizevals=[absdiffsizevals matabsdiffsize(sub2ind(size(matcorrhist), taarrays(:,1), taarrays(:,5+i)))];
	absdifflumvals=[absdifflumvals matabsdifflum(sub2ind(size(matcorrhist), taarrays(:,1), taarrays(:,5+i)))];
	orcorrvals=[orcorrvals matcorror(sub2ind(size(matcorrhist), taarrays(:,1), taarrays(:,5+i)))];
	corrvals1=[corrvals1 matcorr(sub2ind(size(matcorr), taarrays1(:,1), taarrays1(:,5+i)))];
	histcorrvals1=[histcorrvals1 matcorrhist(sub2ind(size(matcorrhist), taarrays1(:,1), taarrays1(:,5+i)))];
	absdiffsizevals1=[absdiffsizevals1 matabsdiffsize(sub2ind(size(matcorrhist), taarrays1(:,1), taarrays1(:,5+i)))];
	absdifflumvals1=[absdifflumvals1 matabsdifflum(sub2ind(size(matcorrhist), taarrays1(:,1), taarrays1(:,5+i)))];
	orcorrvals1=[orcorrvals1 matcorror(sub2ind(size(matcorrhist), taarrays1(:,1), taarrays1(:,5+i)))];
end
absdiffsizevals = absdiffsizevals-repmat(mean(absdiffsizevals')',1,6);[x i] = max(absdiffsizevals'); absdiffsizesel=absdiffsizevals==repmat(x', 1, 6); absdiffsizeobjs=taarrays(:,6:end)'; absdiffsizeobjs=absdiffsizeobjs(:); absdiffsizeobjs=absdiffsizeobjs(absdiffsizesel'==1);
absdifflumvals = absdifflumvals-repmat(mean(absdifflumvals')',1,6);[x i] = max(absdifflumvals'); absdifflumsel=absdifflumvals==repmat(x', 1, 6); absdifflumobjs=taarrays(:,6:end)'; absdifflumobjs=absdifflumobjs(:); absdifflumobjs=absdifflumobjs(absdifflumsel'==1);
orcorrvals = orcorrvals-repmat(mean(orcorrvals')',1,6);[x i] = max(orcorrvals'); orcorrsel=orcorrvals==repmat(x', 1, 6); orcorrobjs=taarrays(:,6:end)'; orcorrobjs=orcorrobjs(:); orcorrobjs=orcorrobjs(orcorrsel'==1);
histcorrvals = histcorrvals-repmat(mean(histcorrvals')',1,6);[x i] = max(histcorrvals'); histcorrsel=histcorrvals==repmat(x', 1, 6); histcorrobjs=taarrays(:,6:end)'; histcorrobjs=histcorrobjs(:); histcorrobjs=histcorrobjs(histcorrsel'==1);
zcorrvals = zscore(corrvals')';
corrvals = corrvals-repmat(mean(corrvals')',1,6);[x i] = max(corrvals'); corrsel=corrvals==repmat(x', 1, 6); corrobjs=taarrays(:,6:end)'; corrobjs=corrobjs(:); corrobjs=corrobjs(corrsel'==1); [i j]=ind2sub(size(corrsel'), find(corrsel'));corrpos=i;
absdiffsizevals1 = absdiffsizevals1-repmat(mean(absdiffsizevals1')',1,6);[x i] = max(absdiffsizevals1'); absdiffsizesel1=absdiffsizevals1==repmat(x', 1, 6); absdiffsizeobjs1=taarrays1(:,6:end)'; absdiffsizeobjs1=absdiffsizeobjs1(:); absdiffsizeobjs1=absdiffsizeobjs1(absdiffsizesel1'==1);
absdifflumvals1 = absdifflumvals1-repmat(mean(absdifflumvals1')',1,6);[x i] = max(absdifflumvals1'); absdifflumsel1=absdifflumvals1==repmat(x', 1, 6); absdifflumobjs1=taarrays1(:,6:end)'; absdifflumobjs1=absdifflumobjs1(:); absdifflumobjs1=absdifflumobjs1(absdifflumsel1'==1);
orcorrvals1 = orcorrvals1-repmat(mean(orcorrvals1')',1,6);[x i] = max(orcorrvals1'); orcorrsel1=orcorrvals1==repmat(x', 1, 6); orcorrobjs1=taarrays(:,6:end)'; orcorrobjs1=orcorrobjs1(:); orcorrobjs1=orcorrobjs1(orcorrsel1'==1);
histcorrvals1 = histcorrvals1-repmat(mean(histcorrvals1')',1,6);[x i] = max(histcorrvals1'); histcorrsel1=histcorrvals1==repmat(x', 1, 6); histcorrobjs1=taarrays1(:,6:end)'; histcorrobjs1=histcorrobjs1(:); histcorrobjs1=histcorrobjs1(histcorrsel1'==1);
corrvals1 = corrvals1-repmat(mean(corrvals1')',1,6);[x i] = max(corrvals1'); corrsel1=corrvals1==repmat(x', 1, 6); corrobjs1=taarrays1(:,6:end)'; corrobjs1=corrobjs1(:); corrobjs1=corrobjs1(corrsel1'==1);[i j]=ind2sub(size(corrsel1'), find(corrsel1'));corrpos1=i;histcorrsel1=double(histcorrsel1); absdffsizesel1=double(absdiffsizesel1); absdifflumsel1=double(absdifflumsel1); orcorrsel1=double(orcorrsel1);corrsel1=double(corrsel1);
corrsel=double(corrsel); modsel=double(modsel);modsel2=double(modsel1);subsel=double(subsel);subsel1=double(subsel1);subsel2=double(subsel2);
%mutualinfo(subpos1, subpos2)
%[x p] = partialcorr(double(subsel(:)), double(modsel(:)), [corrvals(:) orcorrvals(:) absdifflumvals(:)]) % Hmmm, the values are not independent, so the p-value may not mean much
%[x p] = corr(double(subsel(:)), double(modsel(:)))

objsorigorder=[arrays1; arrays1; arrays1; arrays1(allowrep2,:)];
objsorigorder=objsorigorder(sel,:);
taattvalsreord=[];
t=[];
ctrl=[]; ctrl_line=[ 579        1103         820         756         908         493]; % A control for the biases towards certain positions - has no effect, though
for n=1:size(taattvals,1)
	t= [t; circshift(objsorigorder(n,end-5:end), [0 find(taarrays(n,end-5:end) == objsorigorder(n,1))-1])];
	taattvalsreord= [taattvalsreord; circshift(taattvals(n,:), [0 find(taarrays(n,end-5:end) == objsorigorder(n,1))-1])];
	ctrl=[ctrl; ctrl_line];
end
[x i]=sort(taattvalsreord,2);
[x rankattvalsreord]=sort(i,2);
z1=taattvalsreord-repmat(mean(taattvalsreord')', 1, 6);
z2=zscore(taattvalsreord')';
%[x p] = regress(zscore(double(subsel(:))), zscore(double(z2(:))))
%[x p] = regress(zscore(double(subsel(:))), [zscore(double(z2(:))) zscore(double(zcorrvals(:)))])
%[x p] = partialcorr(zscore(double(subsel(:))), zscore(double(z2(:))), zscore(double(zcorrvals(:))))





%We now want to see if the model tends to choose objects that are "similar" to the target, for various measure of "similar".

marraysmod=[tgt1' tgtpos1'  fixedobjsmod];
tp=marraysmod(tgtpos1~=-10,:);
%marrayssubjs= [tgt1' tgtpos1' fixedobjs1];
%tpsub = marrayssubjs(marrayssubjs(:,2)~=-10,:);


distvals=[];histcorrvals=[]; orcorrvals=[]; absdiffsizevals=[]; absdifflumvals=[]; corrvals=[]; distC2vals=[];
for i=1:6
	% ord: sorts the fixated position in order of their pixelwise correlation to the target
	%sim: for the successively fixated objects, compute their average similarity to the target - excluding actual fixations to the target
	distvals=[distvals matdist(sub2ind(size(matdist), tp(:,1), tp(:,2+i)))]; [tmp distord] = sort(distvals,2); simdist{i} = (distvals(distvals(:,i)>0.0,i)); % because 0 means "actual target" 
	distC2vals=[distC2vals matdistC2(sub2ind(size(matdistC2), tp(:,1), tp(:,2+i)))]; [tmp distC2ord] = sort(distC2vals,2); simdistC2{i} = (distC2vals(distC2vals(:,i)>0.0,i)); % because 0 means "actual target" 
	corrvals=[corrvals matcorr(sub2ind(size(matcorr), tp(:,1), tp(:,2+i)))]; [tmp corrord] = sort(corrvals,2); simcorr{i} = (corrvals(corrvals(:,i)<1.0,i)); % because 1.0 means "actual target" 
	histcorrvals=[histcorrvals matcorrhist(sub2ind(size(matcorrhist), tp(:,1), tp(:,2+i)))]; [tmp histcorrord] = sort(histcorrvals,2); simhistcorr{i} =(histcorrvals(histcorrvals(:,i)<1.0,i)); 
	absdiffsizevals=[absdiffsizevals matabsdiffsize(sub2ind(size(matcorrhist), tp(:,1), tp(:,2+i)))]; [tmp absdiffsizeord] = sort(absdiffsizevals, 2); simabsdiffsize{i} = ( absdiffsizevals(absdiffsizevals(:,i)>0,i)); % because 0 = target 
	absdifflumvals=[absdifflumvals matabsdifflum(sub2ind(size(matcorrhist), tp(:,1), tp(:,2+i)))]; [tmp absdifflumord] = sort(absdifflumvals, 2); simabsdifflum{i} = (absdifflumvals(absdifflumvals(:,i)>0.0,i)); 
	orcorrvals=[orcorrvals matcorror(sub2ind(size(matcorrhist), tp(:,1), tp(:,2+i)))]; [tmp orcorrord] = sort(orcorrvals,2); simor{i} =(orcorrvals(orcorrvals(:,i)<1.0,i)); 
end

myf=figure; hold on;
for i=1:6
	subplot(611);hold on; title('Euclidean distance');
	errorbar(i, mean(simdist{i})-2000, std(simdist{i}) / sqrt(length(simdist{i}))); bar(i, mean(simdist{i})-2000);set(gca, 'XTick', [0:7]); set(gca, 'XTickLabel', [' '; '1' ;'2'; '3'; '4'; '5'; '6'; ' ']); 
	subplot(612);hold on; title('Pixelwise correlation');
	errorbar(i, mean(simcorr{i}), std(simcorr{i}) / sqrt(length(simcorr{i}))); bar(i, mean(simcorr{i}));set(gca, 'XTickLabel', [' '; '1' ;'2'; '3'; '4'; '5'; '6'; ' ']);
	subplot(613);hold on; title('Luminance histogram correlation');
	errorbar(i, mean(simhistcorr{i}), std(simhistcorr{i}) / sqrt(length(simhistcorr{i}))); bar(i, mean(simhistcorr{i}));set(gca, 'XTickLabel', [' '; '1' ;'2'; '3'; '4'; '5'; '6'; ' ']);
	subplot(614);hold on; title('Absolute mean difference in luminance'); 
	errorbar(i, mean(simabsdifflum{i}), std(simabsdifflum{i}) / sqrt(length(simabsdifflum{i}))); bar(i, mean(simabsdifflum{i}));set(gca, 'XTickLabel', [' '; '1' ;'2'; '3'; '4'; '5'; '6'; ' ']);
	subplot(615);hold on; title('Absolute difference in size');
	errorbar(i, mean(simabsdiffsize{i}), std(simabsdiffsize{i}) / sqrt(length(simabsdiffsize{i}))); bar(i, mean(simabsdiffsize{i}));set(gca, 'XTickLabel', [' '; '1' ;'2'; '3'; '4'; '5'; '6'; ' ']);
	subplot(616);hold on; title('Distance between C2b vectors');
	errorbar(i, mean(simdistC2{i}), std(simdistC2{i}) / sqrt(length(simdistC2{i}))); bar(i, mean(simdistC2{i})); set(gca, 'XTickLabel', [' '; '1' ;'2'; '3'; '4'; '5'; '6'; ' ']);
end
set(myf, 'Position', [100 100 500 2200]);
 set(gcf,'PaperPositionMode','auto');
%print('-depsc', '-r300', 'modfixdist.eps');
%%print('-dpng', '-r300', 'modfixdist.png');




%vals=[]; 
%for i=1:1000 
%	z=rand(size(taarrays1,1),6); randsel1=double(z==repmat(max(z')', 1, 6));
%	z=rand(size(taarrays1,1),6); randsel2=double(z==repmat(max(z')', 1, 6));
%	vals=[vals  partialcorr(randsel2(:), subsel1(:), [corrsel1(:) orcorrsel1(:) absdifflumsel1(:)])]; 
%end



% Which trials are we going to include in our matrix analyses?
% Note that tgtpos1 == tgtpos2 (normally!)

select1 = tgtpos1' == -10  & fixedobjs1(:,1) ~= -1;
select1rep = tgtpos1rep' == -10  & fixedobjs1rep(:,1) ~= -1;
select2 = tgtpos2' == -10  & fixedobjs2(:,1) ~= -1;
select2rep = tgtpos2rep' == -10  & fixedobjs2rep(:,1) ~= -1;  

%select1 = tgtpos1' ~= -10  & fixedobjs1(:,1) ~= -1 & fixedobjs1(:,1) ~= tgt1' & fixedobjsmod(:,1) ~= tgt1';
%select1rep = tgtpos1rep' ~= -10  & fixedobjs1rep(:,1) ~= -1 & fixedobjs1rep(:,1) ~= tgt1rep'& fixedobjsmod(:,1) ~= tgt1rep';
%select2 = tgtpos2' ~= -10  & fixedobjs2(:,1) ~= -1 & fixedobjs2(:,1) ~= tgt2'& fixedobjsmod(:,1) ~= tgt2';
%select2rep = tgtpos2rep' ~= -10  & fixedobjs2rep(:,1) ~= -1 & fixedobjs2rep(:,1) ~= tgt2rep'& fixedobjsmod(:,1) ~= tgt2rep';

%select1 =  fixedobjs1(:,1) ~= -1;
%select1rep =  fixedobjs1rep(:,1) ~= -1;
%select2 = fixedobjs2(:,1) ~= -1;
%select2rep = fixedobjs2rep(:,1) ~= -1;  

select2rep(1:880) = 0; select2rep(end-1320:end)=0; % VERY IMPORTANT!! Remove the blocks that are just duplicate first sessions due to subjects not taking second session


if EXCLUDE1234 == 1
	select1(1:880)=0; select1rep(1:880)=0; select2(1:880)=0; select2rep(1:880)=0; % If you want to eliminate subjects 1 2 10 & 11 (blocks 1 & 2)
end

%select2rep(1:880) = 0; select2rep(end-440:end)=0; % VERY IMPORTANT!!

%select1 = tgtpos1' ~= -10  & fixedobjs1(:,1) ~= -1;
%select1rep = tgtpos1rep' ~= -10  & fixedobjs1rep(:,1) ~= -1;
%select2 = tgtpos2' ~= -10  & fixedobjs2(:,1) ~= -1;
%select2rep = tgtpos2rep' ~= -10  & fixedobjs2rep(:,1) ~= -1;  
%select2rep(1:880) = 0; select2rep(end-440:end)=0; % VERY IMPORTANT!!

%select1(ceil(end-end/8):end) = 0;
%select2(ceil(end-end/8):end) = 0;

tgnum=[tgt1(select1) tgt1rep(select1rep) tgt2(select2) tgt2rep(select2rep)];
tgpos=[tgtpos1(select1) tgtpos1rep(select1rep) tgtpos2(select2) tgtpos2rep(select2rep)];
arrayssel = [arrays1(select1,:); arrays1rep(select1rep,:); arrays2(select2,:); arrays2rep(select2rep,:)];
attvalssel = [attvals(select1,:); attvals(select1rep,:); attvals(select2,:); attvals(select2rep,:)];
oia=[isobjinarray1(select1, :); isobjinarray1rep(select1rep, :); isobjinarray2(select2, :); isobjinarray2rep(select2rep, :)];
subjs=[subjs1(select1), subjs1rep(select1rep), subjs2(select2), subjs2rep(select2rep)];
subobjs=[fixedobjs1(select1,1); fixedobjs1rep(select1rep,1); fixedobjs2(select2,1); fixedobjs2rep(select2rep,1)]';
modobjs=[fixedobjsmod(select1,1); fixedobjsmod(select1rep,1); fixedobjsmod(select2,1); fixedobjsmod(select2rep,1)]'; %Assuming tgtpos1 = tgtpos2 = tgtposmod, at least for -10!
z= arrayssel  == repmat(modobjs',1,6); [modpos j] = ind2sub( size(z'), find(z'));
z= arrayssel  == repmat(subobjs',1,6); [subpos j] = ind2sub( size(z'), find(z'));

ratio12m=attvalssel(:,1)./attvalssel(:,2); 
errs = subobjs' ~= tgnum'; errm = modobjs' ~= tgnum'; diffms = subobjs' ~= modobjs';
[x p] = corr(errs, ratio12m);



mat=zeros(NBOBJS);
for so=1:NBOBJS
    for to=1:NBOBJS
        if (sum(oia(:, so) & tgnum'==to') > 0)
            mat(so, to) = sum(subobjs==so & tgnum == to) / sum(oia(:, so) & tgnum'==to');
        end;
    end;
end;
matsub=mat;
mat=zeros(NBOBJS);
for so=1:NBOBJS
    for to=1:NBOBJS
        if (sum([isobjinarray1(select1, so); isobjinarray1rep(select1rep,so)] & [tgt1(select1)'; tgt1rep(select1rep)'] ==to) > 0)
            mat(so, to) = sum([fixedobjs1(select1,1); fixedobjs1rep(select1rep,1)] ==so & [tgt1(select1)'; tgt1rep(select1rep)'] == to) / sum([isobjinarray1(select1, so); isobjinarray1rep(select1rep,so)] ...
				& [tgt1(select1)';tgt1rep(select1rep)'] ==to);
        end;
    end;
end;
matsub1=mat;
mat=zeros(NBOBJS);
for so=1:NBOBJS
    for to=1:NBOBJS
        if (sum([isobjinarray2(select2, so); isobjinarray2rep(select2rep,so)] & [tgt2(select2)'; tgt2rep(select2rep)'] ==to) > 0)
            mat(so, to) = sum([fixedobjs2(select2,1); fixedobjs2rep(select2rep,1)] ==so & [tgt2(select2)'; tgt2rep(select2rep)'] == to) / sum([isobjinarray2(select2, so); isobjinarray2rep(select2rep,so)] ...
				& [tgt2(select2)';tgt2rep(select2rep)'] ==to);
        end;
    end;
end;
matsub2=mat;

mat=zeros(NBOBJS);
for so=1:NBOBJS
    for to=1:NBOBJS
        if (sum(oia(:, so) & tgnum'==to') > 0)
            mat(so, to) = sum(modobjs==so & tgnum == to) / sum(oia(:, so) & tgnum'==to');
        end;
    end;
end;
matmod=mat;
mat=zeros(NBOBJS);
for so=1:NBOBJS
    for to=1:NBOBJS
        if (sum([isobjinarray1(select1, so); isobjinarray1rep(select1rep,so)] & [tgt1(select1)'; tgt1rep(select1rep)'] ==to) > 0)
            mat(so, to) = sum([fixedobjsmod(select1,1); fixedobjsmod(select1rep,1)] ==so & [tgt1(select1)'; tgt1rep(select1rep)'] == to) / sum([isobjinarray1(select1, so); isobjinarray1rep(select1rep,so)] ...
				& [tgt1(select1)';tgt1rep(select1rep)'] ==to);
        end;
    end;
end;
matmod1=mat;
mat=zeros(NBOBJS);
for so=1:NBOBJS
    for to=1:NBOBJS
        if (sum([isobjinarray2(select2, so); isobjinarray2rep(select2rep,so)] & [tgt2(select2)'; tgt2rep(select2rep)'] ==to) > 0)
            mat(so, to) = sum([fixedobjsmod(select2,1); fixedobjsmod(select2rep,1)] ==so & [tgt2(select2)'; tgt2rep(select2rep)'] == to) / sum([isobjinarray2(select2, so); isobjinarray2rep(select2rep,so)] ...
				& [tgt2(select2)';tgt2rep(select2rep)'] ==to);
        end;
    end;
end;
matmod2=mat;


% NOTE: this particular section is better done in rgrss.m 
matn=zeros(NBOBJS);matn1=zeros(NBOBJS);matn2=zeros(NBOBJS);
for so=1:NBOBJS
	for to=1:NBOBJS
		matn(so, to) = sum(oia(:, so) & tgnum'==to');
		matn1(so, to) = sum([isobjinarray1(select1, so); isobjinarray1rep(select1rep,so)]  & [tgt1(select1)';tgt1rep(select1rep)'] ==to);	   
		matn2(so, to) = sum([isobjinarray2(select2, so); isobjinarray2rep(select2rep,so)]  & [tgt2(select2)';tgt2rep(select2rep)'] ==to);	   
	end
end;
zn=matn(1:end, 1:end); zn1=matn1(1:end, 1:end);  zn2=matn2(1:end, 1:end); zn12=min(zn1, zn2);
%zn=min(zn1, zn2); % Note: this is necessary to make sure that all have the same number of values, so we can control zs1/2 with the similarity values
zc=matcorr(1:end, 1:end); zc=zc(eye(size(zc)) == 0 & zn>=MINCOUNT); zc=zscore(zc); %We need to take out the diagonal, even if we remove correct trials (because then the lower diagonal induces spuriours corr!) 
zs=matsub(1:end, 1:end); zs=zs(eye(size(zs)) == 0 & zn>=MINCOUNT);  zs=zscore(zs);
zs1=matsub1(1:end, 1:end); zs1=zs1(eye(size(zs1)) == 0 & zn12>=MINCOUNT);  zs1=zscore(zs1);
zs2=matsub2(1:end, 1:end); zs2=zs2(eye(size(zs2)) == 0 & zn12>=MINCOUNT);  zs2=zscore(zs2);
zm=matmod(1:end, 1:end); zm=zm(eye(size(zm)) == 0 & zn>=MINCOUNT); zm=zscore(zm);
zm1=matmod1(1:end, 1:end); zm1=zm1(eye(size(zm1)) == 0 & zn12>=MINCOUNT); zm1=zscore(zm1);
zm2=matmod(1:end, 1:end); zm2=zm2(eye(size(zm2)) == 0 & zn12>=MINCOUNT); zm2=zscore(zm2);
zhist=matcorrhist(1:end, 1:end); zhist=zhist(eye(size(zn)) == 0 & zn>=MINCOUNT);  zhist=zscore(zhist);
zabsdiffsize=matabsdiffsize(1:end, 1:end); zabsdiffsize=zabsdiffsize(eye(size(zn)) == 0 & zn>=MINCOUNT);  zabsdiffsize=zscore(zabsdiffsize);
zabsdifflum=matabsdifflum(1:end, 1:end); zabsdifflum=zabsdifflum(eye(size(zn)) == 0 & zn>=MINCOUNT);  zabsdifflum=zscore(zabsdifflum);
zobjsizes=matobjsizes(1:end, 1:end); zobjsizes=zobjsizes(eye(size(zn)) == 0 & zn>=MINCOUNT); zobjsizes=zscore(zobjsizes);
zobjlums=matobjlums(1:end, 1:end); zobjlums=zobjlums(eye(size(zn)) == 0 & zn>=MINCOUNT); zobjlums=zscore(zobjlums);
zobjcontrasts=matobjcontrasts(1:end, 1:end); zobjcontrasts=zobjcontrasts(eye(size(zn)) == 0 & zn>=MINCOUNT); zobjcontrasts=zscore(zobjcontrasts);
zobjmsdiffbg=matobjmsdiffbg(1:end, 1:end); zobjmsdiffbg=zobjmsdiffbg(eye(size(zn)) == 0 & zn>=MINCOUNT); zobjmsdiffbg=zscore(zobjmsdiffbg);
zobjsizeabsdevs=matobjsizeabsdevs(1:end, 1:end); zobjsizeabsdevs=zobjsizeabsdevs(eye(size(zn)) == 0 & zn>=MINCOUNT); zobjsizeabsdevs=zscore(zobjsizeabsdevs);
zobjlumabsdevs=matobjlumabsdevs(1:end, 1:end); zobjlumabsdevs=zobjlumabsdevs(eye(size(zn)) == 0 & zn>=MINCOUNT); zobjlumabsdevs=zscore(zobjlumabsdevs);
zor=matcorror(1:end, 1:end); zor=zor(eye(size(zn)) == 0 & zn>=MINCOUNT);  zor=zscore(zor);
matdummy=zeros(NBOBJS); matdummy(1:MINCOUNT,:)=1; zdummy=matdummy(1:end, 1:end); zdummy=zdummy(eye(size(zn)) == 0 & zn>=MINCOUNT); zdummy=zscore(zdummy);
 % NOTE: We should NOT include bottom-up factors in this regression!
[b bint] = regress(zs, [ones(size(zs)) zc  zhist zor zabsdiffsize zabsdifflum zm]) 
bmod = b(end);
[b bint] = regress(zs, [ones(size(zs)) zc zm])
[c p] = corr(zs, zm)
[pcorrmodzc p] = partialcorr(zs, zm, zc)
[pcorrmodall p] = partialcorr(zs, zm, [zc zhist zor zabsdiffsize zabsdifflum ])
%[x p] = partialcorr(zs, zm, [zc  zhist zor zobjmsdiffbg zobjcontrasts zobjsizes zabsdiffsize zobjsizeabsdevs zobjlums zobjlumabsdevs]) 
%[b bint] = regress(zs, [ones(size(zs)) zc  zhist zor zobjmsdiffbg zobjcontrasts zobjsizes zabsdiffsize zobjsizeabsdevs zabsdifflum zobjlums zobjlumabsdevs zm]) 



bs=[]; bints=[]; pcorrs=[];
for i=1:1000
	if mod(i,100)==1; disp(i); end
	% We want to generate randomized model responses (i.e. destination of 1st
	% fixation) that nevertheless have the same performance as the original
	% model (i.e. 1st saccade is to the actual target on target-present trials
	% with the same frequency).
	% To do this, we first generate a series of responses that are all wrong,
	% i.e. different from target. Then, on target-present trials, we set the
	% appropriate number of trials to have the correct response
	nbcorrect = sum(fmod(tgtpos1 ~= -10) == tgt1(tgtpos1 ~= -10)'); 
	fwrong=-1*ones(1,size(fixorder,1)); for i=1:size(fixorder,1); fwrong(i)=tgt1(i); while(fwrong(i)==tgt1(i)) uf=ceil(6*rand());  fwrong(i) = arraysmodel(i, uf); end; end; 
	selpos= find(tgtpos1 ~= -10); selpos=selpos(randperm(numel(selpos))); selpos=selpos(1:nbcorrect); 
	frnd=fwrong; frnd(selpos)=tgt1(selpos);
	%fixedobjsmodrnd=-1*ones(size(fixorder)); for i=1:size(fixorder,1); uf=fixorder(ceil(rand()*size(fixorder,1)),:);  if rand()<.6 & tgtpos1(i)~=-10 uf(1)=tgtpos1(i); end; for j=1:numel(uf); if uf(j)==-1 continue; end; fixedobjsmodrnd(i,j) = arraysmodel(i, uf(j));  end; end; % tgtpos1 is supposed to be the same as for the model
	%modobjsrnd=[fixedobjsmodrnd(select1,1); fixedobjsmodrnd(select1rep,1); fixedobjsmodrnd(select2,1); fixedobjsmodrnd(select2rep,1)]'; 
	modobjsrnd=[frnd(select1) frnd(select1rep) frnd(select2) frnd(select2rep)]; 
	mat=zeros(NBOBJS);
	for so=1:NBOBJS
		for to=1:NBOBJS
			if (sum(oia(:, so) & tgnum'==to') > 0)
				mat(so, to) = sum(modobjsrnd==so & tgnum == to) / sum(oia(:, so) & tgnum'==to');
			end;
		end;
	end;
	matmodrnd=mat;
	zmrnd=matmodrnd(1:end, 1:end); zmrnd=zmrnd(eye(size(zn)) == 0 & zn>=MINCOUNT); zmrnd=zscore(zmrnd);
	%[x p] = partialcorr(zs, zmrnd, zc); 
	[x p] = partialcorr(zs, zmrnd, [zc zhist zor zabsdiffsize zabsdifflum ]);
	pcorrs=[pcorrs x];
	%[b bint] = regress(zs, [ones(size(zs)) zc  zhist zor zabsdiffsize zmrnd]); 
	%bs = [bs b(end)]; bints = [bints; bint(end,:)];
end;
%disp(['Proportion of randomized data returning a larger b for model in the multi-factor regression than real data: ' num2str(mean(bs>bmod))]);
disp(['Proportion of randomized data returning a larger partial corr for model than real data: ' num2str(mean(pcorrs>pcorrmodall))]);



% 40x40 subject-model confusion matrix: each x,y cell contains the number of
% trials where subject answered x and model answered y, divided by the total
% number of trials in which both x and y were present
matsm=-1*ones(NBOBJS);
for so=1:NBOBJS
    for mo=1:NBOBJS
        if (sum(oia(:, so) & oia(:,mo)) > 0)
            matsm(so, mo) = sum(subobjs==so & modobjs == mo) / sum(oia(:, so) & oia(:,mo));
        end;
    end;
end;
nodiag= matsm(eye(size(matsm)) == 0);  
disp(['Mean value of subject-model confusion matrix, outside diagonal: ' num2str(mean(nodiag(:))) ', within diagonal: ' num2str(mean(diag(matsm)))]);
disp(['P-value of the anova between these two groups: ' num2str(anova1(matsm(:), reshape(eye(size(matsm)), numel(matsm), 1), 'off'))]);
disp('Note: ONLY USEFUL FOR TARGET-ABSENT TRIALS');


matsubjs={}; matpossubjs={}; matmssubjs={};
for s=1:16
	mat=zeros(NBOBJS);
	for so=1:NBOBJS
		for to=1:NBOBJS
			if (sum(oia(subjs==s, so) & tgnum(subjs==s)'==to') > 0)
				mat(so, to) = sum(subobjs(subjs==s)==so & tgnum(subjs==s) == to) / sum(oia((subjs==s), so) & tgnum(subjs==s)'==to');
			end;
		end;
	end;
	matsubjs{s} = mat;
	mat=zeros(NBOBJS);
	for so=1:NBOBJS
		for mo=1:NBOBJS
			if (sum(oia(subjs==s, so) & oia(subjs==s,mo)) > 0)
				mat(so, mo) = sum(subobjs(subjs==s)==so & modobjs(subjs==s) == mo) / sum(oia(subjs==s, so) & oia(subjs==s,mo));
			end;
		end;
	end;
	matmssubjs{s} = mat;
	mat=zeros(6); 
	for so=1:6
		for mo=1:6
				mat(so, mo) = sum(subpos(subjs==s)==so & modpos(subjs==s)==mo) ;
		end;
	end;
	mat=mat./ repmat(sum(mat')', 1, 6);  
	matpossubjs{s} = mat;
end


matpos=zeros(6); 
for i=1:numel(tgnum); matpos(subpos(i), modpos(i)) =  matpos(subpos(i), modpos(i)) +1; end;
matpos=matpos./ repmat(sum(matpos')', 1, 6);  













if 1 == 0
	figure;mat=zeros(6); e=eye(6);
	for s=1:16
		subplot(2,8,s);
		imagesc(matpossubjs{s});
		disp(num2str(anova1(matpossubjs{s}(:), e(:), 'off')));
		%imagesc(matsubjs{s});
		%imagesc(matmssubjs{s});
		mat=mat+matpossubjs{s};
	end
end





return;









% The rest of this code is unused.


%trialst=1;trialend=1+400;
%firsts=firsts(trialst:trialend);tgtpos=tgtpos(trialst:trialend);fixorder=fixorder(trialst:trialend,:);usattvals=usattvals(trialst:trialend,:);uscorrwtgt=uscorrwtgt(trialst:trialend,:); fixedpos=fixedpos(trialst:trialend,:);
%isobjinarray=isobjinarray(trialst:trialend, :); tgt=tgt(trialst:trialend); arrays=arrays(trialst:trialend,:); usdists = usdists(trialst:trialend);attvals=attvals(trialst:trialend,:);

% On Target-absent trials, first subj == first mod much more often than chance
mat=[firsts' fixorder(:,1) tgtpos'];
matta = mat(mat(:,end)==-10 , :);
mattp = mat(mat(:,end)~=-10 , :); matbe = mat(mat(:,end)~=-10 & mat(:,1) ~= mat(:,end) & mat(:,2) ~= mat(:,end) , :);
%modsub = mean(matbe(:,1) == matbe(:,2));
modsub = mean(matta(:,1) == matta(:,2));
for n=1:10000
	matr=[firsts(randperm(2000))' fixorder(randperm(2000),1) tgtpos(randperm(2000))'];
	%matber = matr(matr(:,end)~=-10 & matr(:,1) ~= matr(:,end) & matr(:,2) ~= matr(:,end) , :);
	mattar = matr(matr(:,end)==-10 , :);
	modsub = [modsub ; mean(mattar(:,2) == mattar(:,1))];
end
disp([modsub(1) mean(modsub>modsub(1))]);


% Human and model performance
% Success rates after each fixation (excluding target absent trials)
fo=fixorder(tgtpos~=-10, :);
so=fixedpos(tgtpos~=-10, :);
tgp=tgtpos(tgtpos~=-10);
cumsum(mean(so==repmat(tgp',1,6)))
cumsum(mean(fo==repmat(tgp',1,6)))
% Average number of fixations




% On ALL TRIALS, correlation between attention value and being subject first choice (excluding actual target objects) - using Euclidean distance or correlation for similarity 

f=firsts;
tg=tgtpos;
a=usattvals;
%a = a - repmat(min(a')', 1, 6);
%a = a ./ repmat(max(a')', 1, 6);
%a=zscore(a')';
%asim=usdists;
asim=uscorrwtgt;
%asim = asim - repmat(min(asim')', 1, 6);
%asim = asim ./ repmat(max(asim')', 1, 6);
%asim=zscore(asim')';
isf=zeros(size(a)); 
fidx=sub2ind(size(a), 1:size(a,1), f); isf(fidx)=1;
istg=zeros(size(a));
% We want to set to 1 those positions were the actual target was, using tg to tell us which column it is for each line - but not at these lines where tg==-10!
tsubs=[ 1:size(a,1); tg]; tsubs=tsubs(:, tg~=-10);
tidx=sub2ind(size(a), tsubs(1,:), tsubs(2,:)); istg(tidx)=1;
anottg=a(istg==0);
asimnottg=asim(istg==0);
isfnottg=isf(istg==0);
%kruskalwallis(asimnottg, isfnottg, 'off')
[x p] = corr(anottg(:), isfnottg(:), 'type', 'spearman')
[x p] = partialcorr(anottg(:), isfnottg(:), asimnottg(:), 'type', 'spearman')
%errorbar([0 1], [mean(anottg(isfnottg==0)) mean(anottg(isfnottg==1))], [std(anottg(isfnottg==0))/sqrt(numel(anottg(isfnottg==0))) std(anottg(isfnottg==1))/sqrt(numel(anottg(isfnottg==1)))])



% In error and target-absent trials, correlation between attention value and being subject first choice (excluding actual target objects MIGHT BE BIASED) - using Euclidean distant for similarity 

firstsr=firsts;%(randperm(numel(firsts)));
f=firstsr(firstsr~=tgtpos);
tg=tgtpos(firstsr~=tgtpos);
a=usattvals(firstsr~=tgtpos, :);
%a = a - repmat(min(a')', 1, 6);
%a = a ./ repmat(max(a')', 1, 6);
a=zscore(a')';
asim=usdists(firstsr~=tgtpos, :);
%asim=uscorrwtgt(firstsr~=tgtpos, :);
%asim = asim - repmat(min(asim')', 1, 6);
%asim = asim ./ repmat(max(asim')', 1, 6);
%asim=zscore(asim')';
isf=zeros(size(a)); 
fidx=sub2ind(size(a), 1:size(a,1), f); isf(fidx)=1;
istg=zeros(size(a));
% We want to set to 1 those positions were the actual target was, using tg to tell us which column it is for each line - but not at these lines where tg==-10!
tsubs=[ 1:size(a,1); tg]; tsubs=tsubs(:, tg~=-10);
tidx=sub2ind(size(a), tsubs(1,:), tsubs(2,:)); istg(tidx)=1;
anottg=a(istg==0);
asimnottg=asim(istg==0);
isfnottg=isf(istg==0);
%kruskalwallis(asimnottg, isfnottg, 'off')
[x p] = corr(anottg(:), isfnottg(:), 'type', 'spearman')
[x p] = partialcorr(anottg(:), isfnottg(:), asimnottg(:), 'type', 'spearman')







% Target-absent trials
f=firsts(tgtpos==-10);
tg=tgtpos(tgtpos==-10);
a=usattvals(tgtpos==-10, :);
a=zscore(a')';
%asim=usdists(tgtpos==-10, :);
asim=uscorrwtgt(tgtpos==-10, :);
asim=zscore(asim')';
isf=zeros(size(a)); 
fidx=sub2ind(size(a), 1:size(a,1), f); isf(fidx)=1;
[x p] = partialcorr(a(:), isf(:), asim(:))




% All trials where both subject 1st and model 1st differ from target - including errors by both and target-absent trials
% Probably most reliable (but less data)
f=firsts(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)');
tg=tgtpos(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)');
a=usattvals(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)',:);
a=zscore(a')';
asim=uscorrwtgt(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)', :);
%asim=usdists(tgtpos~=firsts & tgtpos~=fixorder(:,1)', :);
asim=zscore(asim')';
isf=zeros(size(a)); 
fidx=sub2ind(size(a), 1:size(a,1), f); isf(fidx)=1;
[x p] = corr(a(:), isf(:))
[x p] = partialcorr(a(:), isf(:), asim(:))



% Correlation between rank given by model and being subjet first choice (you want that to be negative if you don't zero out!)
f=firsts(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)');
tg=tgtpos(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)');
fo=fixorder(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)', :);
%fo=fo(randperm(size(fo,1)),:);
[asim fosim] = sort(uscorrwtgt(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)', :), 2, 'descend');
co = uscorrwtgt(firsts ~= 5 & tgtpos~=firsts & tgtpos~=fixorder(:,1)', :);
%f=firsts(tgtpos~=firsts & tgtpos~=fixorder(:,1)');
%tg=tgtpos(tgtpos~=firsts & tgtpos~=fixorder(:,1)');
%fo=fixorder(tgtpos~=firsts & tgtpos~=fixorder(:,1)', :);
%[asim fosim] = sort(uscorrwtgt(tgtpos~=firsts & tgtpos~=fixorder(:,1)', :), 2, 'descend');
[x rankmodel]= sort(fo, 2, 'ascend');
[x ranksim]= sort(fosim, 2, 'ascend');
%asim=usdists(tgtpos~=firsts & tgtpos~=fixorder(:,1)', :);
rankmodel(rankmodel~=1)=0;
%ranksim(ranksim~=2)=0; % Really 2! You know 1st ranksim won't be 1st choice, since 1st ranksim is always the target! But wait, not true for the target-absent trials... But zeroing out ranksim eliminates a lot of information from it! 
isf=zeros(size(fo)); 
fidx=sub2ind(size(fo), 1:size(fo,1), f); isf(fidx)=1;
[x p] = corr(rankmodel(:), isf(:))
[x p] = corr(ranksim(:), isf(:))
[x p] = partialcorr(rankmodel(:), isf(:), ranksim(:))
%[x p] = corr(co(:), isf(:))
%[x p] = partialcorr(rankmodel(:), isf(:), co(:))




%
ar=arrays(tgtpos~=firsts,: );
a=attvals(tgtpos~=firsts,  :);
fo=fixorder(tgtpos~=firsts, :);
oia=isobjinarray(firsts ~=tgtpos, :);
f=firsts(firsts ~= tgtpos );
tgp=tgtpos(firsts~=tgtpos );
tgnum=tgt(firsts ~=tgtpos );
[asim fosim] = sort(uscorrwtgt(firsts ~=tgtpos,:), 2, 'descend');

ar=arrays(tgtpos~=firsts & fixorder(:,1)' ~= tgtpos, :);
a=attvals(tgtpos~=firsts & fixorder(:,1)' ~= tgtpos, :);
fo=fixorder(tgtpos~=firsts & fixorder(:,1)' ~= tgtpos, :);
oia=isobjinarray(firsts ~=tgtpos & fixorder(:,1)' ~= tgtpos, :);
f=firsts(firsts ~= tgtpos & fixorder(:,1)' ~= tgtpos);
tgp=tgtpos(firsts~=tgtpos & fixorder(:,1)' ~= tgtpos);
tgnum=tgt(firsts ~=tgtpos & fixorder(:,1)' ~= tgtpos);
[asim fosim] = sort(uscorrwtgt(firsts ~=tgtpos & fixorder(:,1)' ~= tgtpos,:), 2, 'descend');


ar=arrays;
a=attvals;
fo=fixorder;
oia=isobjinarray;
f=firsts;
tgp=tgtpos;
tgnum=tgt;
[asim fosim] = sort(uscorrwtgt, 2, 'descend');

ar=arrays(tgtpos == -10,: );
a=attvals(tgtpos == -10,  :);
fo=fixorder( tgtpos == -10, :);
oia=isobjinarray( tgtpos == -10, :);
f=firsts( tgtpos == -10);
tgp=tgtpos( tgtpos == -10);
tgnum=tgt(tgtpos == -10);
[asim fosim] = sort(uscorrwtgt(tgtpos == -10,:), 2, 'descend');



ar=arrays(tgtpos ~= -10,: );
a=attvals(tgtpos ~= -10,  :);
fo=fixorder( tgtpos ~= -10, :);
oia=isobjinarray( tgtpos ~= -10, :);
f=firsts( tgtpos ~= -10);
tgp=tgtpos( tgtpos ~= -10);
tgnum=tgt(tgtpos ~= -10);
[asim fosim] = sort(uscorrwtgt(tgtpos ~= -10,:), 2, 'descend');


ar=arrays(firsts~=5 & tgtpos ~= -10,: );
a=attvals(firsts~=5 & tgtpos ~= -10,  :);
fo=fixorder(firsts~=5 & tgtpos ~= -10, :);
oia=isobjinarray(firsts ~= 5 & tgtpos ~= -10, :);
f=firsts(firsts ~= 5 & tgtpos ~= -10);
tgp=tgtpos(firsts ~= 5  & tgtpos ~= -10);
tgnum=tgt(firsts ~= 5 & tgtpos ~= -10);
[asim fosim] = sort(uscorrwtgt(firsts ~=5 & tgtpos ~= -10,:), 2, 'descend');

ar=arrays(firsts~=5,: );
a=attvals(firsts~=5,  :);
fo=fixorder(firsts~=5, :);
oia=isobjinarray(firsts ~= 5, :);
f=firsts(firsts ~= 5);
tgp=tgtpos(firsts ~= 5 );
tgnum=tgt(firsts ~= 5);
[asim fosim] = sort(uscorrwtgt(firsts ~=5,:), 2, 'descend');

ar=arrays(tgtpos~=firsts & fixorder(:,1)' ~= tgtpos & firsts~=5, :);
a=attvals(tgtpos~=firsts & fixorder(:,1)' ~= tgtpos & firsts~=5, :);
fo=fixorder(tgtpos~=firsts & fixorder(:,1)' ~= tgtpos & firsts~=5, :);
oia=isobjinarray(firsts ~=tgtpos & fixorder(:,1)' ~= tgtpos & firsts~=5, :);
f=firsts(firsts ~= tgtpos & fixorder(:,1)' ~= tgtpos & firsts~=5);
tgp=tgtpos(firsts~=tgtpos & fixorder(:,1)' ~= tgtpos & firsts~=5);
tgnum=tgt(firsts ~=tgtpos & fixorder(:,1)' ~= tgtpos & firsts~=5);
[asim fosim] = sort(uscorrwtgt(firsts ~=tgtpos & fixorder(:,1)' ~= tgtpos & firsts~=5 ,:), 2, 'descend');



%resp=[]; resx=[];for n=1:1000 % This randomization does not apply to the partial correlation though, since it breaks any correlation between firsts and similarity... But it allows us to see if there is a bias in the correlation p-value computation?
%firstsr=firsts; %if (n>1); firstsr=firsts(randperm(numel(firsts))); end;
% ar=arrays(tgtpos~=firstsr & fixorder(:,1)' ~= tgtpos, :);
% a=attvals(tgtpos~=firstsr & fixorder(:,1)' ~= tgtpos, :);
% fo=fixorder(tgtpos~=firstsr & fixorder(:,1)' ~= tgtpos, :);
% oia=isobjinarray(firstsr ~=tgtpos & fixorder(:,1)' ~= tgtpos, :);
% f=firstsr(firstsr ~= tgtpos & fixorder(:,1)' ~= tgtpos);
% tgp=tgtpos(firstsr~=tgtpos & fixorder(:,1)' ~= tgtpos);
% tgnum=tgt(firstsr ~=tgtpos & fixorder(:,1)' ~= tgtpos);
% [asim fosim] = sort(uscorrwtgt(firstsr ~=tgtpos & fixorder(:,1)' ~= tgtpos,:), 2, 'descend');
% Target - subject response confusion matrix, divided at each point by the number of trials where the target object was target and the responded object was present
f2 = f; %f2=f2(randperm(numel(f2)));
fidx=sub2ind(size(ar), 1:size(ar,1), f2); subjobjs=ar(fidx);
mat=zeros(NBOBJS);
for so=1:NBOBJS
	for to=1:NBOBJS
		if (sum(oia(:, so) & tgnum'==to') > 0)
			mat(so, to) = sum(subjobjs==so & tgnum == to) / sum(oia(:, so) & tgnum'==to');
		end;
	end;
end;
matsub=mat;
% Target - model response confusion matrix, divided at each point by the number of trials where the target object was target and the responded object was present
f2 = fo(:,1)'; %f2=f2(randperm(numel(f2)));
fidx=sub2ind(size(ar), 1:size(ar,1), f2); modobjs=ar(fidx);
mat=zeros(NBOBJS);
for so=1:NBOBJS
	for to=1:NBOBJS
		if (sum(oia(:, so) & tgnum'==to') > 0)
			mat(so, to) = sum(modobjs==so & tgnum == to) / sum(oia(:, so) & tgnum'==to');
		end;
	end;
end;
matmod=mat;
% Similarity-based response confusion matrix, divided at each point by the number of trials where the target object was target and the responded object was present
f2 = fosim(:,1)'; %f2=f2(randperm(numel(f2)));
fidx=sub2ind(size(ar), 1:size(ar,1), f2); modobjs=ar(fidx);
mat=zeros(NBOBJS);
for so=1:NBOBJS
	for to=1:NBOBJS
		if (sum(oia(:, so) & tgnum'==to') > 0)
			mat(so, to) = sum(modobjs==so & tgnum == to) / sum(oia(:, so) & tgnum'==to');
		end;
	end;
end;
matsim=mat;
matn=zeros(NBOBJS);
for so=1:NBOBJS
	for to=1:NBOBJS
		matn(so, to) = sum(oia(:, so) & tgnum'==to');
	end;
end;
zn=matn(1:end, 1:end); 
zc=matcorr(1:end, 1:end); zc=zc(eye(size(zc)) == 0 & zn>=2); %We need to take out the diagonal, even if we remove correct trials (because then the lower diagonal induces spuriours corr!) 
zs=matsub(1:end, 1:end); zs=zs(eye(size(zs)) == 0 & zn>=2);  
zm=matmod(1:end, 1:end); zm=zm(eye(size(zm)) == 0 & zn>=2);
zsim=matsim(1:end, 1:end); zsim=zsim(eye(size(zsim)) == 0 & zn>=2);
zhist=matcorrhist(1:end, 1:end); zhist=zhist(eye(NBOBJS) == 0 & zn>=2); 
zabsdiffsize=matabsdiffsize(1:end, 1:end); zabsdiffsize=zabsdiffsize(eye(NBOBJS) == 0 & zn>=2); 
zobjsizes=matobjsizes(1:end, 1:end); zobjsizes=zobjsizes(eye(NBOBJS) == 0 & zn>=2); 
zobjcontrasts=matobjcontrasts(1:end, 1:end); zobjcontrasts=zobjcontrasts(eye(NBOBJS) == 0 & zn>=2); 
zobjmeanabs=matobjmeanabs(1:end, 1:end); zobjmeanabs=zobjmeanabs(eye(NBOBJS) == 0 & zn>=2); 
zor=matcorror(1:end, 1:end); zor=zor(eye(NBOBJS) == 0 & zn>=2); 
matdummy=zeros(NBOBJS); matdummy(1:2,:)=1; zdummy=matdummy(1:end, 1:end); zdummy=zdummy(eye(NBOBJS) == 0 & zn>=2); 
[b bint] = regress(zs, [ones(size(zs)) zc  zhist zor zobjcontrasts./1000 zobjsizes./1000 zm])  % NOTE: CONFOUND! SOME OF THESE BOTTOM-UP FACTORS ARE SIZE-DEPENDENT! 

 [x p]= corr(zs, zm, 'type', 'spearman') ;
 [x p]= partialcorr(zs, zm, zsim, 'type', 'spearman') 
 [x p]= partialcorr(zs, zm, zc, 'type', 'spearman')   
 %  resx=[resx; x]; resp=[resp; p];
% disp(n); end;



 % 40x40 subject-model confusion matrix: each x,y cell contains the number of
% trials where subject answered x and model answered y, divided by the total
% number of trials in which both x and y were present
%f=f(randperm(numel(f)));
mat=zeros(NBOBJS);
matc=zeros(NBOBJS);
for i=1:length(f)
	mat(ar(i,f(i)), ar(i,fo(i,1))) = mat(ar(i,f(i)), ar(i,fo(i,1))) + 1;
end;
for subjc=1:NBOBJS
	for modc=1:NBOBJS
		mat(subjc, modc) = mat(subjc, modc) / sum(oia(:, subjc) & oia(:, modc));
		matc(subjc, modc) = sum(oia(:, subjc) & oia(:, modc));
	end
end
	nodiag= mat(matc > 1 & eye(size(mat)) == 0);  s = mean(nodiag(:))
mean(diag(mat))


res=[];
for n=1:1000
	f2=f(randperm(numel(f))); % If we use trials with both-errors, this is a bias, because both-errors exclude the target from possible choices for each trial and this randomization doesn't
	matr=zeros(NBOBJS);
	matc=zeros(NBOBJS);
	for i=1:length(f2)
		matr(ar(i,f2(i)), ar(i,fo(i,1))) = matr(ar(i,f2(i)), ar(i,fo(i,1))) + 1;
	end;
	for subjc=1:NBOBJS
		for modc=1:NBOBJS
			matr(subjc, modc) = matr(subjc, modc) / sum(oia(:, subjc) & oia(:, modc));
			matrc(subjc, modc) = sum(oia(:, subjc) & oia(:, modc));
		end
	end
	%nodiag= matr(matrc > 1 & eye(size(matr)) == 0);  res=[res; mean(nodiag(:))];
	res=[res; mean(diag(matr))];
end;

%e=eye(size(mat)); e=e(:); kruskalwallis(mat(:), e, 'off') % kruskal-wallis test seems to give 'significant' results even for randomized firsts...




%% How easy is each object to find?
%easysubj=[]; easymod=[];
%for i=1:NBOBJS
%	easysubj=[easysubj sum(tgt==i & tgtpos ~= -10 & tgtpos==firsts)/sum(tgt==i & tgtpos ~= -10)];
%	easymod=[eagymod sum(tgt==i & tgtpos ~= -10 & tgtpos==fixorder(:,1)')/sum(tgt==i & tgtpos ~= -10)];
%end;
%diso=[]; 
%fgr i=1:NBOBJS
%		diso = [diso; (load(['./out/C2bout_ndp_obj' num2str(i) '.png.raw.out']))'];
%end;
%


% Position matrix. Not biased, but must correct for the subject's different base probabilities of choosing different positions:
z=zeros(6); for i=1:numel(f); z(f(i), fo(i,1)) =  z(f(i), fo(i,1)) +1; end;
z=z./ repmat(sum(z')', 1, 6);  
imagesc(z);





