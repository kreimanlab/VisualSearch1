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




