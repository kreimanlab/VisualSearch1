function []=runModelArgs(varargin)

NBOBJS = 40; % some problem with object 40, img 12 or 13?...
NBPICSPEROBJ = 40;  % multiple of 4!
NBNAT=500;
UNITTYPE = 'ndp';
NBPROTS = 600;
SIGMA_IMG = 600;
SIGMA_ATT = 5;

RANDOMIZE=0;
IMAGETYPE = 'x';
for myarg=varargin
    if strcmp(myarg, 'array') == 1
        IMAGETYPE = 'ARRAY';
    end
    if strcmp(myarg, 'unpred') == 1
        IMAGETYPE = 'UNPRED';
    end
    if strcmp(myarg, 'rand') == 1
        RANDOMIZE = 1;
    end
end

if length(IMAGETYPE) < 3
    error('Must provide at least one argument "array" or "unpred"')
end

tic;
disp(';;;;');

img={};imgn={};imgnofb={};
img1={}; img2={}; 
img1n={}; img2n={}; 
img1nofb={}; img2nofb={}; 
img_g={};imgnofb_g={};

if strcmp(IMAGETYPE, 'QUADSIZE')
	uposfile = load('./images/quadsize.txt');
elseif strcmp(IMAGETYPE, 'ARRAY')
	uposfile = load('./images/array.txt');
elseif strcmp(IMAGETYPE, 'QUAD')
	uposfile = load('./images/quad.txt');
elseif strcmp(IMAGETYPE, 'UNPRED')
	uposfile = load('./images/unpredpos.txt');
else 
	error('Which image type?');
end;
ux={}; uy={}; usize={};
for i=1:size(uposfile,1)
	ux{uposfile(i,1)}{uposfile(i,2)} = uposfile(i, 3);
	uy{uposfile(i,1)}{uposfile(i,2)} = uposfile(i, 4);
	usize{uposfile(i,1)}{uposfile(i,2)} = uposfile(i, 5);
end;

%d=[]; 
%for i=1:NBOBJS
%	disp(['Cluttered data for object ' num2str(i)]);
%	for j=1:NBPICSPEROBJ
%		if strcmp(IMAGETYPE, 'QUADSIZE')
%			d = [d; (load(['./out/C2bout_ndp_' num2str(i) '.quadsize' num2str(j) '.ot.png.raw.out']))'];
%		elseif strcmp(IMAGETYPE, 'UNPRED')
%			d = [d; (load(['./out/C2bout_ndp_' num2str(i) '.unpred' num2str(j) '.ot.png.raw.out']))'];
%		else 
%			error('Which image type?');
%		end;
%	end;
%end;

diso=[];
for i=1:NBOBJS
		disp(['Isolated image for object ' num2str(i)]);
		diso = [diso; (load(['./out/C2bout_ndp_' num2str(i) '.small.ot.png.raw.out']))'];
end;

dnat=[];

%for i=NBNAT/2:NBNAT   
for i=1:NBNAT  
		disp(['Data for natural image ' num2str(i)]);
%		dnat = [dnat; (load(['./out/C2bout_ndp_nat' num2str(i) '.png.raw.out']))'];
end;

coeffs={}; origdiffmeans={};



for NUMOBJ=1:NBOBJS % 1:NBOBJS

 % Computing the top-down modulation coefficient for each object

		m1=diso(NUMOBJ, :);
%		m2=mean(dnat);
		m2=mean(diso);
		

		diffmeans = m1 ./  m2; %- m2;
		
		origdiffmeans{NUMOBJ}=diffmeans;
		

		
		diffmeans = diffmeans - min(diffmeans(:)); %zscore(diffmeans); 
		diffmeans(diffmeans < 0 ) = 0; %s=sort(diffmeans); diffmeans(diffmeans < s(floor(end-length(diffmeans)/20)))=0;
		diffmeans = diffmeans ./  max(diffmeans(:)); 
		diffmeans = 1 +  diffmeans;


		coeffs{NUMOBJ}  = 1 * diffmeans; 


end;



for NUMOBJ=1:NBOBJS 

	for numim=1:NBPICSPEROBJ
		
		mycoeffs = coeffs{NUMOBJ};
		
		othero=ceil(rand()*NBOBJS);
		while (othero == NUMOBJ) othero = ceil(rand()*NBOBJS); end;
                
                if RANDOMIZE > 0
		    mycoeffs = (coeffs{othero}); % Chooses top-down parameters from a random object, for each input image, as a randmized control
                end
		
		disp(['Loading S2b outputs for image ' num2str(numim) ' of object ' num2str(NUMOBJ)]);
		if strcmp(IMAGETYPE, 'QUADSIZE')
				fid = fopen(['./out/S2bout_' UNITTYPE '_' num2str(NUMOBJ) '.quadsize' num2str(numim) '.ot.png.raw.out.short'], 'r');
		elseif strcmp(IMAGETYPE, 'ARRAY')
				%fid = fopen(['/groups/kreiman/thomas/out/S2bout_' UNITTYPE '_' num2str(NUMOBJ) '.array' num2str(numim) '.ot.png.raw.out.short'], 'r');
				fid = fopen(['./out/S2bout_' UNITTYPE '_' num2str(NUMOBJ) '.array' num2str(numim) '.ot.png.raw.out.short'], 'r');
		elseif strcmp(IMAGETYPE, 'QUAD')
				fid = fopen(['./out/S2bout_' UNITTYPE '_' num2str(NUMOBJ) '.quad' num2str(numim) '.ot.png.raw.out.short'], 'r');
		elseif strcmp(IMAGETYPE, 'UNPRED')
				fid = fopen(['./out/S2bout_' UNITTYPE '_' num2str(NUMOBJ) '.unpred' num2str(numim) '.ot.png.raw.out.short'], 'r');
		else 
			error('Which image type?');
		end;
		data = fread(fid, 'uint16');
		fclose(fid);
		
		p1 = data(1:5:end) + 1;  % +1 because Matlab array indices range from 1 to MAX, while C ranges from 0 to MAX-1
		x1 = data(4:5:end) ;
		y1 = data(3:5:end) ;
		v1 =  double(data(5:5:end));
		v1 = v1 ./ 65535.0;
	
	%sub - sub: Found 0.37375 times (chance level .25). dif/dif : same

		%v1 = (v1 - sumS);% ./ stdS;
		%v1 = (v1 ./ (sumS+.01)); %./ stdS;

		v1(v1 < 0) = 0;
		v1 = v1 + 1e-15;
		
		thisimgfb=zeros(256); thisimgnofb=zeros(256); 
		for j=1:numel(p1) % For some reason, putting numel(v1) makes it 10 times as slow!!!
			thisimgnofb(x1(j), y1(j)) = thisimgnofb(x1(j), y1(j)) + v1(j);   
			thisimgfb(x1(j), y1(j)) = thisimgfb(x1(j), y1(j)) + v1(j) * mycoeffs(p1(j));  % +1;
		end;
		
		
		sumfbg = zeros(256); sumnofbg = zeros(256);

		img{NUMOBJ}{numim} = thisimgfb;
		imgnofb{NUMOBJ}{numim} = thisimgnofb + 1e-15;

		img_n = thisimgfb ./ (thisimgnofb + SIGMA_IMG);
		imgnofb_n = thisimgnofb ./ (thisimgnofb + SIGMA_IMG);

		imgn{NUMOBJ}{numim} = img{NUMOBJ}{numim} ./ imgnofb{NUMOBJ}{numim};
		%imgsn{NUMOBJ}{numim} = img_n ./ (imgnofb_n + SIGMA_ATT) ; 
		imgsn{NUMOBJ}{numim} = img{NUMOBJ}{numim} ./ (imgnofb{NUMOBJ}{numim} + SIGMA_ATT) ;  
	end;
	disp('Done');

end;

nbcells=zeros(256)+1e-15;
for j=1:numel(p1)
	nbcells(x1(j), y1(j)) = nbcells(x1(j), y1(j)) + 1;
end;
nbcells(71,71)=1e-15;nbcells(111,71)=1e-15;nbcells(151,71)=1e-15; nbcells(71,111)=1e-15;nbcells(111,111)=1e-15;nbcells(151,111)=1e-15; nbcells(71,151)=1e-15;nbcells(111,151)=1e-15;nbcells(151,151)=1e-15;
smoothnbcells = imfilter(nbcells, fspecial('gaussian', 121, 42), 'replicate') ;

%NBOBJS=20; NBPICSPEROBJ=20;
IORRAD = 50;  ERRMARG=0; MAXTRIES=5; IORFILTER = fspecial('gaussian', IORRAD*2+1, IORRAD/3);IORFILTER = 1 - .2*(IORFILTER ./ max(IORFILTER(:)));
FALLRATERAD=25;
found=[]; maxqs = []; maxcs=[]; maxls=[]; maxvs=[];quadrs=[];nbtries=[]; frates=[]; rsizes=[];  pics={}; oamaps={}; amaps={}; numamap=1;
nbtriesperobj = [];
for o=1:NBOBJS
	nbtriesthisobj = [];
	for q=1:NBPICSPEROBJ
		if strcmp(IMAGETYPE, 'QUADSIZE'); pic = imread(['./images/' num2str(o) '.quadsize' num2str(q) '.ot.png']);	
		elseif strcmp(IMAGETYPE, 'QUAD'); pic = imread(['./images/' num2str(o) '.quad' num2str(q) '.ot.png']);
		elseif strcmp(IMAGETYPE, 'ARRAY'); pic = imread(['./images/' num2str(o) '.array' num2str(q) '.ot.png']);
		elseif strcmp(IMAGETYPE, 'UNPRED'); pic = imread(['./images/' num2str(o) '.unpred' num2str(q) '.ot.png']);
		else error('Which image type?');
		end;
		amap = imgsn{o}{q} ;
		%amap = imgsn{ceil(rand()*NBOBJS)}{ceil(rand()*NBPICSPEROBJ)} ;
		amap(71,71)=1e-15;amap(111,71)=1e-15;amap(151,71)=1e-15; amap(71,111)=1e-15;amap(111,111)=1e-15;amap(151,111)=1e-15; amap(71,151)=1e-15;amap(111,151)=1e-15;amap(151,151)=1e-15;
		%amap=imfilter(amap, fspecial('gaussian', 15, 2), 'replicate') ./ imfilter(nbcells, fspecial('gaussian', 15, 2), 'replicate') ; t=zeros(size(amap));t(8:246,8:246)=amap(8:246,8:246);amap=t; amap2=amap; 
		z=zeros(256); 
		% TO CALCULATE FRATE: 
		%z=imfilter(amap, fspecial('gaussian', 121, 42), 'replicate') ./ smoothnbcells ; t=zeros(size(z));t(8:246,8:246)=z(8:246,8:246);z=t;  % larger filter = better perf!
		%s=sort(amap(:)); amap(amap<s(end-200)) = s(end-200);
		
		amaps{numamap} = amap;  pics{numamap} = pic;
		fallrate{o}{q}=1e-15;
		hasfound=0;
		[maxv maxn] = max(amap(:));[myl myc] = ind2sub(size(amap), maxn); maxcs=[maxcs; myc]; maxls=[maxls; myl];
		for numtry=1:MAXTRIES

			[maxv maxn] = max(amap(:));[myl myc] = ind2sub(size(amap), maxn); 
			if (myc+ERRMARG > ux{o}{q} & myc-ERRMARG < ux{o}{q} + usize{o}{q} - 1 ...
				& myl+ERRMARG > uy{o}{q} & myl-ERRMARG < uy{o}{q} + usize{o}{q} - 1)
				hasfound = 1;
				outer = z(max(1,myl-FALLRATERAD):min(256, myl+FALLRATERAD), max(1,myc-FALLRATERAD):min(256, myc+FALLRATERAD));
				inner = z(max(1,myl-(FALLRATERAD-1)):min(256, myl+(FALLRATERAD-1)), max(1,myc-(FALLRATERAD-1)):min(256, myc+(FALLRATERAD-1)));
				fallrate{o}{q} = maxv / ((sum(outer(:))-sum(inner(:)))/(length(outer(:))-length(inner(:)))+1e-15);
				break;
			else
				t=zeros(256+2*IORRAD); t(IORRAD+1:end-IORRAD, IORRAD+1:end-IORRAD) = amap;
				t(myl:myl+2*IORRAD, myc:myc+2*IORRAD) =  t(myl:myl+2*IORRAD, myc:myc+2*IORRAD) .* IORFILTER;
				amap = t(IORRAD+1:end-IORRAD, IORRAD+1:end-IORRAD) ;
			end;
		end;
		
		if hasfound == 0 % if target not found within maximum allowed number of tries (=fixations), number of tries is recorded as -1
			numtry = -1;
		end;
		nbtries=[nbtries numtry];
		nbtriesthisobj=[nbtriesthisobj numtry];
		found = [found hasfound];
		
		frates=[frates; fallrate{o}{q}];
		rsizes=[rsizes; usize{o}{q}];
		numamap = numamap+1;
	end;
	nbtriesperobj = [nbtriesperobj; nbtriesthisobj]; 
end;
foundattry=zeros(1,MAXTRIES); foundattryperobj=zeros(NBOBJS,MAXTRIES); foundthroughtry=zeros(1,MAXTRIES); foundthroughtryperobj=zeros(NBOBJS,MAXTRIES); 
for t=1:MAXTRIES
	foundattry(t) = (length(nbtries(nbtries==t)) / length(nbtries)); 
	foundthroughtry(t) = sum(foundattry(1:t));
	for o=1:NBOBJS
		tmp = nbtriesperobj(o,:);
		foundattryperobj(o,t) = (length(tmp(tmp==t)) / length(tmp)); % length(tmp) should be = NBPICSPEROBJ
		foundthroughtryperobj(o, t) = sum(foundattryperobj(o, 1:t));
	end;
end;
%disp(['Found ' num2str(foundattry(1)) ' times 1st, ' num2str(sum(foundattry(1:2))) ' times 2nd, ' num2str(sum(foundattry(1:3))) ' times 3rd, ' num2str(sum(foundattry(1:4))) ' times 4th, '   num2str(sum(foundattry(1:5))) ' times 5th.']); 
%disp(['Found ' num2str(foundthroughtry(1)) ' times 1st, ' num2str(foundthroughtry(2)) ' times 2nd, ' num2str(foundthroughtry(3)) ' times 3rd, ' num2str(foundthroughtry(4)) ' times 4th, '   num2str(foundthroughtry(5)) ' times 5th.']); 
disp(['Found ' num2str(mean(foundthroughtryperobj(:,1))) ' times 1st, ' num2str(mean(foundthroughtryperobj(:,2))) ' times 2nd, ' num2str(mean(foundthroughtryperobj(:,3))) ' times 3rd, ' num2str(mean(foundthroughtryperobj(:,4))) ' times 4th, '    num2str(mean(foundthroughtryperobj(:,5))) ' times 5th.']); 
disp(['SEM: ' num2str(std(foundthroughtryperobj(:,1)) / sqrt(NBOBJS)) ' times 1st, ' num2str(std(foundthroughtryperobj(:,2)) / sqrt(NBOBJS)) ' times 2nd, ' num2str(std(foundthroughtryperobj(:,3)) / sqrt(NBOBJS)) ' times 3rd, ' num2str(std(foundthroughtryperobj(:,4)) / sqrt(NBOBJS)) ' times 4th, '   num2str(std(foundthroughtryperobj(:,5)) / sqrt(NBOBJS)) ' times 5th.']); 
fr=frates(found==1); rs=rsizes(found==1); [x, p] = (corr(fr, rs, 'type' , 'pearson'));

return;

%toc;
%set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'Bold');
%z=nbtries; z(z==-1)=6; x=[]; for i=1:NBOBJS; x = [x sum(z(1+(i-1)*NBPICSPEROBJ:i*NBPICSPEROBJ))'./NBPICSPEROBJ]  ; end; figure; 
%plot(x', mean(diso')','*'); ylabel('Average feature activity (isolated image)'); xlabel('Average number of fixations (max. 6)'); axis([1 6 .4 .71]);[r p] = corr(x',mean(diso')') 
%[b bint]=regress(mean(diso')', [x' ones(size(x'))]); hold on; plot(1:.01:6, b(2)+b(1)*(1:.01:6), 'r');
%% To try with and without normalization, (same pattern for normal and random weights): r = -.073, p= .654 / r = -.8731, p < 2*10e-13  --     text(4.5,.68,'\fontsize{16}r = -0.8731');text(4.5,.66,'\fontsize{16}p < 2*10e-13');
%title('With Normalization');
%print('-depsc', '-r300', 'corrDisoFixNorm.eps');


%
%for i=1:20; disp(mean(found(1+(i-1)*40:i*40))); end;  % Useful to see which objects are identified more or less often. Should compare with non-rand coefficients.
%for i=1:20; disp([mean(nbtries(1+(i-1)*40:i*40)) std(nbtries(1+(i-1)*40:i*40))] ); end;  % Useful to see which objects are identified more or less often. Should compare with non-rand coefficients.
%x=[]; for i=1:NBOBJS; x = [x (found(1+(i-1)*NBPICSPEROBJ:i*NBPICSPEROBJ))']  ; end;  % Useful to see which objects are identified more or less often. Should compare with non-rand coefficients.
%x=[];for i=1:20; x=[x; ([mean(found(1+(i-1)*40:i*40)) median(origdiffmeans{i})] )]; end;

%%1st try, 20 objs, 20 pics:
%numpic=5; z=amaps{numpic}; z=imfilter(z,fspecial('gaussian',30,4))./imfilter(nbcells,fspecial('gaussian',30,4)); z(z<median(z(:)))=median(z(:));  figure; imagesc(z);figure; imagesc(pics{numpic});colormap(gray); 
%numpic=311; z=amaps{numpic}; z=imfilter(z,fspecial('gaussian',30,4))./imfilter(nbcells,fspecial('gaussian',30,4)); z(z<median(z(:)))=median(z(:));  figure; imagesc(z);figure; imagesc(pics{numpic});colormap(gray); 



