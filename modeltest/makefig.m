
array=[ 0.525 , 0.73 , 0.855 , 0.92125 , 0.95312 ];
arraysem=[0.050072 0.042254 0.035998 0.029439 0.025448 ];
arrayr=[0.10125 0.21812 0.32188 0.42 0.52437]; 
arrayrsem=[ 0.015835 0.024871 0.031435 0.036053 0.038434];

arraynonorm=[ 0.11313  0.235 0.3325  0.42937  0.495 ];
arraynonormsem=[ 0.036361 0.052439 0.060521 0.067123 0.068285 ];

%using normalization parameter 'LocalMax':
%salarray=[0.12063 0.2275 0.3125 0.37688 0.44063];
%salarraysem=[0.038029 0.053938 0.059707 0.064682 0.067146];

%using normalization parameter 'None':
salarray=[0.11938 0.23187 0.33125 0.43937 0.51625];
salarraysem=[0.037334 0.052372 0.059139 0.066064 0.068801];



%unpred=[0.39875   0.5425 0.63312 0.71125  0.775];
%unpredr=[0.14813 0.26625 0.38125 0.46437 0.545];
%salunpred=[0.19875 , 0.3075 , 0.38562 , 0.4525 , 0.505];
unpred=[0.39875 0.54313 0.635 0.71562 0.7775];
unpredsem=[0.041341 0.039609 0.036469 0.031796 0.027074]; 
unpredr=[0.15437 0.26562 0.36312 0.45562 0.52688];
unpredrsem=[0.011494 0.018007 0.01957 0.020617 0.021379];
salunpred=[0.125 0.20938 0.28438 0.35313 0.42];
salunpredsem=[0.01854 0.026322 0.030562 0.033219 0.034763];


unprednonorm=[ 0.17313 0.26625 0.3475 0.42688 0.47688 ];
unprednonormsem=[ 0.032674 0.037867 0.03995 0.041258 0.041161 ];



quadsize=[0.7425 0.92 0.9675 0.99875 0.99875];
quadsizer=[0.26375 0.48 0.6975 0.94875 0.98875];
salquadsize=[0.2575, 0.4925 , 0.65625 , 0.75438 , 0.80063];

%(note: Saliency map for unpred, randomly shuffling the saliency maps between objects: Found 0.0825 times 1st, 0.1475 times 2nd, 0.19563 times 3rd, 0.245 times 4th, 0.29563 times 5th.)
%(note: Saliency map for quadsize, randomly shuffling the saliency maps between objects: Found 0.23625 times 1st, 0.48062 times 2nd, 0.64875 times 3rd, 0.73875 times 4th, 0.78875 times 5th.)
%(note: Saliency map for array, randomly shuffling the saliency maps between objects: Found 0.12125 times 1st, 0.22812 times 2nd, 0.31937 times 3rd, 0.3925 times 4th, 0.4425 times 5th.);


f1=figure;
set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'Bold');
errorbar( array, arraysem, 'ko-', 'linewidth', 1);
hold on; errorbar(arraynonorm, arraynonormsem, 'ko--', 'linewidth', .5);
hold on; errorbar(arrayr, arrayrsem, 'ko:', 'linewidth', .5);
hold on; errorbar( salarray, salarraysem, 'k*:', 'linewidth', .5, 'markersize', 10);
axis([.9 5.1 0 1.1]);
set(gca,'XTick',(1:5));
xlabel('\fontsize{16}Attentional fixations');
ylabel('\fontsize{16}Success rate over all images');
title('\fontsize{16}Composite images - 9 objects');

print( f1, '-r300', '-dpng', 'ResultsArrayPrint')
print( f1, '-r300', '-depsc', 'ResultsArrayPrint.eps')


f2=figure;
set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontWeight', 'Bold');
errorbar(unpred, unpredsem, 'ko-', 'linewidth', 1);
hold on; errorbar(unprednonorm, unprednonormsem, 'ko--', 'linewidth', .5);
hold on; errorbar(unpredr, unpredrsem, 'ko:', 'linewidth', .5);
hold on; errorbar(salunpred, salunpredsem, 'k*:', 'linewidth', .5, 'markersize', 10);
%legend('Normal', 'Randomized weights');
legend('Normal', 'No normalization', 'Randomized', 'Saliency model', 'location', 'NorthWest');
axis([.9 5.1 0 1.1]);
set(gca,'XTick',(1:5));
xlabel('\fontsize{16}Attentional fixations');
ylabel('\fontsize{16}Success rate over all images');
title('\fontsize{16}Natural backgrounds');

print( f2, '-r300', '-dpng', 'ResultsNaturalBackgroundsPrint')
print( f2, '-r300', '-depsc', 'ResultsNaturalBackgroundsPrint.eps')


return;

f3=figure;
plot(quadsize, 'ko-');
hold on; plot(quadsizer, 'ko--');
hold on; plot(salquadsize, 'k+:');
Legend('Normal', 'Randomized weights', 'Saliency', 'Location', 'SouthEast');
set(gca,'XTick',(1:5));
xlabel('\fontsize{16}Attentional fixations');
ylabel('\fontsize{16}Success rate over all images');
title('\fontsize{16}Composite images - 4 objects, varying size');
set(gca, 'FontSize', 16);



