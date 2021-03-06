clc
clear
close all

input = [0
8.72144e-005
0.00130822
0.00252922
0.00392465
0.00505843
0.00750044
0.00470958
0.003663
0.00287807
0.00270365
0.0210187
0.260684
0.119309
0.0457003
0.0432583
0.0537241
0.0532008
0.0558172
0.0475318
0.0411652
0.0368045
0.0130822
0.0113379
0.0162219
0.0300017
0.0428223
0.0448282
0.0225013
0.00645386
0.000174429
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
];


database_data = [0
0.000336134
0.0010084
0.00142857
0.00302521
0.00529412
0.00579832
0.0129412
0.0120168
0.00109244
8.40336e-005
0.0534454
0.221681
0.13521
0.0429412
0.0605042
0.0632773
0.0445378
0.0372269
0.0647899
0.0357983
0.0442857
0.0133613
0.00983193
0.0155462
0.0141176
0.0254622
0.0610924
0.0138655
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
];


figure,
bar(input);set(gca, 'XLim', [0.0 35.0]);set(gca, 'YLim', [0.0 0.265]);%view([90,-30,60]);set(gcf,'renderer','zbuffer'); %set(gca, 'ZLim', [0.0 600000.0]);
set (gcf,'Position',[400,200,350,262], 'color','w')
%title('Distance Histogram Descriptors of The Rotated Model'); 
xlabel('Distance'); ylabel('Frequency');
figure,
bar(database_data);set(gca, 'XLim', [0.0 35.0]);set(gca, 'YLim', [0.0 0.265]);%view([90,-30,60]);set(gcf,'renderer','zbuffer');  %set(gca, 'ZLim', [0.0 600000.0]);
set (gcf,'Position',[400,200,350,262], 'color','w')
%title('Distance Histogram Descriptors of The Origin Model'); 
xlabel('Distance'); ylabel('Frequency');

