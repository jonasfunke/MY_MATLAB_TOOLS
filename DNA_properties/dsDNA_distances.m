function [ d_out ] = dsDNA_distances( n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    
% c5-c5-distances from model.it, straight B-form DNA, H28-50
d = [      1.0000   16.1662
    2.0000   17.6341
    3.0000   18.4325
    4.0000   19.0702
    5.0000   20.2286
    6.0000   22.4586
    7.0000   25.7857
    8.0000   29.8575
    9.0000   34.1909
   10.0000   38.3616
   11.0000   42.0696
   12.0000   45.2423
   13.0000   48.0025
   14.0000   50.5630
   15.0000   53.2050
   16.0000   56.1589
   17.0000   59.5341
   18.0000   63.2463
   19.0000   67.1143
   20.0000   70.9614
   21.0000   74.5725
   22.0000   77.8951
   23.0000   80.9562
   24.0000   83.8806
   25.0000   86.8160
   26.0000   89.9364
   27.0000   93.3185
   28.0000   96.9155
   29.0000  100.6293
   30.0000  104.3212
   31.0000  107.8932
   32.0000  111.2470
   33.0000  114.4048
   34.0000  117.4604
   35.0000  120.5327
   36.0000  123.7377
   37.0000  127.1112
   38.0000  130.6522
   39.0000  134.2823
   40.0000  137.9003
   41.0000  141.4348
   42.0000  144.8007
   43.0000  148.0197
   44.0000  151.1438
   45.0000  154.2905
   46.0000  157.5203
   47.0000  160.9070
   48.0000  164.4150
   49.0000  167.9868
   50.0000  171.5661];

%{
    % distances from 3D-Dart
    d = [1	16.867
2	18.434
3	19.09
4	19.318
5	19.888
6	21.549
7	24.571
8	28.623
9	33.11
10	37.487
11	41.392
12	44.681
13	47.422
14	49.854
15	52.302
16	55.066
17	58.304
18	61.99
19	65.936
20	69.888
21	73.622
22	77.014
23	80.071
24	82.92
25	85.754
26	88.759
27	92.05
28	95.626
29	99.383
30	103.16
31	106.8
32	110.2
33	113.37
34	116.38
35	119.37
36	122.48
37	125.8
38	129.32
39	132.98
40	136.66
41	140.24
42	143.65
43	146.87
44	149.97
45	153.05
46	156.22
47	159.55
48	163.04
49	166.64
50	170.27];
%}
if max(n)>50
    p = polyfit(d(40:50,1), d(40:50,2), 1); %fit last 10 residuals
    n_add = 51:ceil(max(n));
    d = [d ; n_add' p(1).*n_add'+p(2)];
end

d_out = spline(d(:,1), d(:,2), n);


    %plot(1:50, dsDNA_distances(1:50), '+', 0:100,dsDNA_distances(0:100), '.' )

end

