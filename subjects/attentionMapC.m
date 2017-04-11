NBOBJS = 40;
NBNAT=500;
UNITTYPE = 'ndp';
SIGMA_ATT = 5  

PSIZE = 0;

datadir = './subjectsdata/';
results=[]; sequence=[];
z=load([datadir 'ResultsSaccades_AndrewT_blk1.mat']); % Note: the choice of which to use must be synchronized with makepics.m
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_AlissaT_blk2.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_Fetemeh_blk3.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_Mario_blk4.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_Vivian_blk5.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_Sam_blk6.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_Claudio_blk7.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_Alexandru_blk8.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_Mona_blk9.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];
z=load([datadir 'ResultsSaccades_HC_blk10.mat']);
[x i]=sort(z.sequence); results=[results z.results(i)]; sequence=[sequence x];


params = z.params;
BLOCKSIZE = numel(z.sequence);

tgt=[]; for i=1:length(results);  tgt=[tgt results{i}.target];  end
arrays=[]; for i=1:length(results);  arrays=[arrays; results{i}.probes];  end


params.ARRAYSIZE=6;

load('destrect.mat');
assert(length(destrect) == params.ARRAYSIZE);
assert(PSIZE > 0);

% SHOULD WE ZSCORE?...

diso=[]; imgdata={};
for i=1:NBOBJS
	imgdata{i}=imread(['./images/obj' num2str(i) '.png']);
disp(['Isolated image for object ' num2str(i)]);
diso = [diso; (load(['./out/C2bout_ndp_obj' num2str(i) '.png.raw.out']))'];
end;

%dnat=[];
%
%for i=1:NBNAT
%disp(['Data for natural image ' num2str(i)]);
%dnat = [dnat; (load(['./out/C2bout_ndp_img_' num2str(i) '.pgm.raw.out']))'];
%end;

coeffs={};

%t=diso(randperm(size(diso,1)),:); diso=t;  % RANDOMIZATION

for NUMOBJ=1:NBOBJS

		m1=diso(NUMOBJ, :);
		%m2=mean(dnat) + 1e-14;
		m2=mean(diso) + 1e-14;
		%m2=mean(diso([1:NUMOBJ-1 NUMOBJ+1:end],:)) + 1e-14; % Not really different
		
		
		diffmeans = m1 ./  m2; %- m2;

		diffmeans = diffmeans - min(diffmeans(:)); %zscore(diffmeans); 
		diffmeans(diffmeans < 0 ) = 0; %s=sort(diffmeans); diffmeans(diffmeans < s(floor(end-length(diffmeans)/20)))=0;
		diffmeans = diffmeans ./  max(diffmeans(:)); 
		diffmeans = 1 +  diffmeans;



		coeffs{NUMOBJ}  = 1 * diffmeans; 

end;


usattvals=[]; 
fixorder=[];
uscorrwtgt=[]; 
fixordersim=[];
attvals=[]; 
corrwtgt=[];
usdists=[];


for numtrial=1:length(results)
	
		TARGET = tgt(numtrial);

%	TARGET = ceil(rand()*params.ARRAYSIZE);


	mycoeffs = coeffs{TARGET};
	
	othero=ceil(rand()*NBOBJS);
	while (othero == TARGET) othero = ceil(rand()*NBOBJS); end;
	%mycoeffs = (coeffs{othero}); % ACTUAL GOOD RANDOMIZATION - others are not randomizing enough and may present wild fluctuations between runs.
	
	disp(['Loading S2b outputs of trial ' num2str(numtrial)]);
	fid = fopen(['./out/S2bout_ndp_arraySub_blk' num2str(1+floor((numtrial-1) / BLOCKSIZE)) '_num' num2str(sequence(numtrial)) '.png.raw.out.short'], 'r');

		data = fread(fid, 'uint16');
		fclose(fid);
		
		p1a = data(1:5:end) + 1;  % +1 because Matlab array indices range from 1 to MAX, while C ranges from 0 to MAX-1
		s1a = data(2:5:end) + 1;
		x1a = data(4:5:end);
		y1a = data(3:5:end);
		v1a =  double(data(5:5:end));
		v1a = v1a ./ 65535.0;

		p1=p1a;
		v1=v1a;
		x1=x1a;
		y1=y1a;


		v1(v1 < 0) = 0;
		v1 = v1 + 1e-15;
		
		thisimgfb=zeros(PSIZE); thisimgnofb=zeros(PSIZE); 
		nbcells=zeros(PSIZE);
		for j=1:numel(p1) % For some reason, putting numel(v1) makes it 10 times as slow!!!
			thisimgnofb(x1(j), y1(j)) = thisimgnofb(x1(j), y1(j)) + v1(j);   
			thisimgfb(x1(j), y1(j)) = thisimgfb(x1(j), y1(j)) + v1(j) * mycoeffs(p1(j));  % +1;
			nbcells(x1(j), y1(j)) = nbcells(x1(j), y1(j)) + 1;
		end;
		
		thisimgfb(nbcells>median(unique(nbcells))) = 0;  % Usually only 3 values for each prototype and position - 0, 1 or 2
		thisimgnofb(nbcells>median(unique(nbcells))) = 0;
		
		img = thisimgfb;
		imgnofb = thisimgnofb + 1e-15;

		mapsn = img ./ (imgnofb + SIGMA_ATT) ;  

	maxattobj = zeros(1, params.ARRAYSIZE);
	for p=1:params.ARRAYSIZE
		z = mapsn(destrect{p}(1):destrect{p}(2), destrect{p}(3):destrect{p}(4));
		maxattobj(p) = max(z(:));
	end

	atts = maxattobj;
	usattvals=[usattvals; atts]; 
	
	dists=[];
	for numo=1:params.ARRAYSIZE
			 dists = [dists sqrt(sum((double(imgdata{results{numtrial}.probes(numo)}(:)) - double(imgdata{TARGET}(:))).^2))];
	end;
	correls=[];
	for numo=1:params.ARRAYSIZE
			 correls = [correls corr(double(imgdata{results{numtrial}.probes(numo)}(:)), double(imgdata{TARGET}(:)))];
	end;

	attssim = correls;
	uscorrwtgt=[uscorrwtgt; attssim]; 
	usdists = [usdists; dists];

end;

	[attvals fixorder]= sort(usattvals, 2, 'descend')


	tgtmodel = tgt; arraysmodel=arrays;
	%save(['modeloutCdiam9size44_1200_20s.mat'], 'usattvals',  'uscorrwtgt', 'fixorder', 'attvals', 'usdists', 'tgtmodel', 'arraysmodel');
	save(['modelout.mat'], 'usattvals',  'uscorrwtgt', 'fixorder', 'attvals', 'usdists', 'tgtmodel', 'arraysmodel');
return;







