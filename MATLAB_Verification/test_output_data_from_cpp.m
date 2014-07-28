clc
clear
close all

origin = [0.0795775 0.238732 0.397887 0.557042 0.716197 0.875352 1.03451 1.19366 1.35282 1.51197 1.67113 1.83028 1.98944 2.14859 2.30775 2.4669 2.62606 2.78521 2.94437 3.10352 3.26268 3.42183 3.58099 3.74014 3.8993 4.05845 4.21761 4.37676 4.53592 4.69507 4.85423 6.09959e-027 
9.62887 16.8761 12.5767 12.8382 13.3012 9.74014 8.70536 9.27201 8.94741 9.95113 11.9252 13.2306 13.7106 15.3581 19.6887 25.227 32.1437 36.1351 31.4834 30.9325 40.6667 40.5296 32.973 39.3182 45.9567 39.4533 36.4853 37.5455 35.5153 42.7694 58.7813 0.0021196 
45.8366 54.0355 45.8064 48.5175 32.9342 39.3914 46.514 40.201 29.2461 32.5056 26.2855 32.1211 27.9378 52.3556 68.8877 58.0247 46.8143 47.1326 51.1289 67.2921 106.835 74.3914 86.7526 98.8134 85.5928 67.2468 77.9561 135.422 101.67 122.563 129.836 0.985862 
548.209 355.443 293.682 218.357 357.839 615.851 392.132 373.243 234.422 278.283 317.757 428.795 280.392 428.129 294.613 310.024 352.237 254.688 380.322 402.881 344.361 346.102 278.195 263.098 220.01 262.19 414.592 478.992 401.963 416.911 463.916 17.6973 
4136.76 2974.47 1846.14 2315.89 576.098 1749.33 1770.59 1832.54 1198.79 1183.81 2226.61 2102.16 1618.69 1374.09 991.147 1380.61 1319.66 1562.03 1782.82 1395.87 1062.65 1361.2 1392.53 1298.62 1469.98 1510.92 1364.71 1415.6 1091.34 1162.32 1469.46 61.4668 
6925.23 4754.46 714.563 2557.95 1226.27 1080.46 597.976 683.752 2176.13 2810.62 2128.76 1767.65 2315.16 3631.53 3266.8 2047.06 964.177 1618.99 1180.76 1778.9 1846.9 1872.9 2248.64 2242.44 1968.71 1472.37 2028.52 2265.18 2133.64 2147.53 1565.67 23.9801 
10313.2 7760.2 504.353 2599.32 2135.16 1460.36 1141.53 1981.7 1117.46 2289.13 2581.27 2190.09 1289.68 1322.65 1609.48 1940.77 2599.39 3311.04 2777.64 1978.95 1615.07 1931.4 1577.95 1183.78 1154.11 1632.83 2010.45 1841.04 2642.33 2265.42 1868.61 24.5809 
28457.2 20577.2 5725.12 7228.06 4256.35 4156.97 2863.66 5297.46 2298.53 2466.18 4316.81 7993.5 6312.19 10227.9 6226.49 6130.19 5136.84 4967.29 4570.59 5315.57 5974.18 5649.42 4559.09 3242.27 4469.1 4513.95 4510.87 5372 6817.51 5795.81 5926.76 65.2708 
50929.6 39273.9 17997.9 18198 5955.42 10390.7 6683 10322.2 12850 9775.39 10885.6 19205.4 12268.4 15421.1 13807.4 14150.8 12231.2 11820 8113.53 14186.8 13980 11563.3 12621.6 11844.7 9553.87 10783.7 12327.3 10135.6 14874.7 12980.6 9711.2 1032.09 
76270.3 60050.7 21224.6 34313.2 14154.9 18907.3 18629.1 17735.8 12098.1 15165.7 13566.9 11656.9 14491.4 14253 11405.6 21321.9 19279.2 13651 15563.7 15338.8 14725 21671 19855.2 15799.2 17371.7 16611.8 15244.4 14727.5 12620.7 7886.52 8718.84 836.912 
115357 112026 14123.9 26975.6 27761.9 23213.7 25552.3 26918.4 38323.7 53039.9 37321.1 26824.3 28829.4 28839.5 25184.7 28117.3 25428.9 19523.1 20900.7 26699.5 26245.2 23253.9 23600.6 20757.4 15213.6 15583.6 22114.2 25018.7 24789.8 24837.4 25016.8 636.353 
112500 120746 19116.1 12455.1 14399.4 17855.6 20748.9 17424.1 27559.8 47888.9 41062.6 23178 17140.6 16426.1 17824.1 21217.7 16860.3 19295.9 21934.3 20286.5 19921.6 20613.8 18475.1 18368.5 16318.7 18150.4 17622.3 14382.4 15146.5 20198.8 23785.8 1172.38 
104328 90560.6 13469.3 12279.5 14144.3 18919.4 11871.5 28163.4 9883.01 24014.4 19230.5 19874.1 40019 22635.8 26972 20977.2 11711.7 14859.6 15458.7 19796.4 12999.2 22660.8 19895.5 18985.4 25442.8 16731.9 14423.5 13810.8 13367.8 11574.3 13534.2 770.623 
150451 120011 10987.5 33928.2 5488.77 25492.4 41683.5 38350.6 23276.6 19141.5 24439.9 43244.5 39926.5 31790.5 46673.3 30297.4 31287.5 47246.4 31735.1 30614.7 26509.8 25742.1 41650.7 44401.2 37572 31869.2 28961.7 38539.4 45765.5 43797.7 32624.6 189.578 
172662 158887 24458.3 41064.1 19473.7 13728 36209.9 71191.4 47220.4 39430.7 32528.9 31772 26041.3 26652.1 53364.2 41336.2 25137.5 20688.1 34272.8 36260.1 18104.8 26270.1 40189 24281.7 19369.4 29348.4 28249.2 24111.9 21685.7 32794.4 26820.9 550.442 
97870.8 106468 19794.4 29296.5 26887 18593.6 22097.3 25527.1 26862.9 22737.4 15420.2 16322.1 18651.6 12169.9 8795.02 13229 21177.1 27506.9 22872.1 19660.1 31136.4 35626.5 26388.4 21644.4 29427.4 31420.7 17000.4 12083 20668.8 20908.6 13637.8 1025.4 
106895 114485 35423 43282.1 42740.1 47376 43886.1 36918.6 35635.7 36595.2 27952 29894.4 32914.8 35999.4 39828.4 37858.5 39497.6 43710.1 42274 34543.5 29506.7 37533.4 43548.7 33476.4 21329.6 22976.8 30024 29527.4 23670.9 22452 29196.1 140.287 
110054 92905.7 59839.5 75236.8 43469.2 36347.4 18013.2 33939.8 31255.7 25345.7 20940.4 32524.4 26578.7 26847.7 30457.9 29407.2 30100.1 24692.3 35423.9 32266.8 21207.4 26767.3 20844.6 30144 31266 36397.9 36809.2 27387.1 27159 32486.8 28484 18.2456 
155972 172517 126969 108909 44499.2 39172.1 39432.3 58742.5 42818.4 34477.2 28418.7 26233.6 22792.4 18260.7 23482.3 31675.4 29773.9 31339.4 24619.6 22780.9 34819.4 59232.2 55638.5 42390.1 27854.5 32696 39306 34298.3 29410.5 29308.7 40310.3 4539.89 
156864 153142 133658 113028 56283 47028.4 51075.7 45380.8 30790.3 32179.5 46094.9 78340.4 28879.5 32517.1 20014.1 28763.7 32777.3 39421.6 25352 37187.7 46543.2 66923.2 73955 84838.2 84032.1 75615 77180.5 74129.3 59776.4 70584.7 69858.2 1350.39 
152427 132041 146596 149293 82573.2 49160.7 39304.9 37337.5 42118 43692.3 53062 116652 49387.3 43926.4 35006.1 52118.7 39268 41233.9 47051.2 45975.2 37216.8 46607.8 61620.1 62888.6 67285.8 77361.8 64861.6 49155.9 44067.4 54918 52649.4 2106.57 
114019 90743.2 117082 129071 60818.9 33873.7 17398.1 21837.2 12917 19369.3 25449.1 58965.4 30443.9 29878.2 18441.7 30045.9 23318.7 24046.8 25974.4 25224.7 18473.2 21765.7 26611.8 29813.7 33089.4 40536 38252.5 24374 21653.2 31009.6 45787.1 2053.24 
99110.2 61163.5 131122 113771 75230.3 36705.2 47011.7 46859.6 33875.2 43689.5 35547.3 55216.1 35408.6 45982.3 23200.2 47280.6 24787.8 53179.2 33076.3 49896.1 38213.7 45874.8 48137.8 54613.8 45987.1 54012.1 57645.7 46634 48655.9 42280 41521.4 2976.44 
87901.4 65168.3 131726 115893 68451.9 35742 38745.2 46960.8 26257.5 28650.8 36571.6 38766.3 40397.1 32951.5 36870.2 32594.6 25895.3 34480.1 31761.7 44372.7 42684.6 41741.7 33819.3 38521.9 27186 38647 29658.9 36996.7 31470.5 49136.4 33490.6 3634.56 
81821.2 47826.4 153299 119026 75909.2 67819.1 41114.3 91591.8 30342.5 56367.9 33357.6 42581.8 35896.1 37290.5 30676.4 28833.6 24947.4 25685.3 24383 29150.8 33963.1 33362.5 44693 32705.4 39709.3 38456.8 41057.1 47744.9 46737.9 64313.8 47410.3 5150.31 
27233.4 18501.6 45883.7 54367.5 15684.4 17698.8 10221.4 24197.3 12394.7 10757.5 21800.7 15166.4 14917.8 15368.7 20030.5 16567.8 8952 13921.6 14907.5 11029.2 14238.9 16841.2 14461.7 12706.4 15901.8 16438.1 13416.3 15649.8 19267.2 20245.7 20254.1 2100.65 
25673.6 18237 45903.8 46429.3 20800.4 16066 19708.8 35469 16972.9 22777.6 21517.8 15609.8 16712.8 18667.7 21937.9 20252.3 19872.3 19184.1 22485.9 27874.3 26668.3 26605.6 25266.5 24605 21276.5 18354.3 15479.9 15875.6 15002 12960.6 17845.4 6003.72 
10894.2 5671.94 24839.5 19575 10600 11508.4 7735.99 9391.11 5105.96 8005.23 5116.43 4168.89 3523.59 3545.06 4039.14 3881.24 5887.57 5311.06 5480.16 10084.1 3820.63 5326.79 5346.52 3449.04 2222.26 4318.66 4864.52 3589.82 3599.34 5828.58 5919.36 349.078 
4854.94 3220.59 11110.5 13898 6870.78 10562.4 7075.28 8002.59 4214.74 5324.74 8137.51 1664.98 5872.22 5158.29 2604.23 4974.16 3376.65 4000.95 3840.98 3198.56 3169.34 2988.69 2468.02 1676.41 2615.74 1908.39 1705.41 2024.65 1805.83 1800.95 1425.95 222.909 
998.22 2441.29 3893.85 4735.06 4469.35 4428.52 3624.22 2727.52 2279.39 1587.53 1673.92 1645.17 1705.69 2315.32 2111.92 2261.56 2404.61 1875.77 1944.32 1715.08 1415.13 1532.68 1279.07 1272.2 1232.07 974.453 951.546 719.196 564.979 523.172 418.378 36.5403 
779.939 2260.85 3516.39 4432.72 4943.95 5039.32 4761.3 4195.41 3454.01 2657.29 1915.02 1312.04 899.493 692.635 674.594 804.583 1028.4 1289.05 1535.68 1729.68 1847.53 1880.79 1833.81 1720.27 1559.31 1371.94 1178.27 995.373 836.115 708.493 615.592 103.458 
911.082 2670.78 4249.61 5547.96 6494.78 7053.76 7225.22 7043.56 6570.91 5888.06 5083.94 4245.16 3447.05 2747.03 2181.14 1763.62 1489.27 1337.89 1279.9 1282.06 1312.67 1345.4 1361.68 1351.36 1312.09 1247.73 1166.17 1077.1 990.141 913.364 852.459 121.66 
];

change30 = [0.0795775 0.238732 0.397887 0.557042 0.716197 0.875352 1.03451 1.19366 1.35282 1.51197 1.67113 1.83028 1.98944 2.14859 2.30775 2.4669 2.62606 2.78521 2.94437 3.10352 3.26268 3.42183 3.58099 3.74014 3.8993 4.05845 4.21761 4.37676 4.53592 4.69507 4.85423 6.09959e-027 
9.62887 15.4648 10.3534 13.108 15.1279 10.2933 9.29377 12.2629 12.7418 10.8067 9.23272 9.59961 10.4508 12.345 19.0654 29.3037 35.6907 32.1055 28.0161 37.1639 41.5412 28.6278 33.8079 58.5882 53.7866 31.5327 34.728 40.6289 30.1716 32.9986 54.306 0.0207388 
45.8366 73.8058 49.0126 43.2146 36.8959 33.477 36.7587 35.8686 33.8698 27.3142 26.1449 31.5966 38.087 45.2491 49.9044 56.9734 53.1263 47.7266 72.1235 97.5197 92.3582 70.5091 60.9231 73.8574 76.3497 69.4278 90.0013 118.443 125.313 127.482 132.345 2.53941 
548.209 708.979 255.905 168.675 305.855 499.757 406.935 320.394 243.732 204.004 281.235 312.633 408.243 363.476 270.834 270.839 314.645 398.302 422.489 391.92 306.937 325.169 344.564 274.768 266.955 303.217 402.915 415.648 392.044 455.256 539.265 13.3601 
4136.76 5906.84 2083.33 509.505 603.518 1372.6 1869.45 1522.82 1098.44 1387.37 2075.25 2010.41 1271.34 930.716 1083.79 1390.28 1740.41 1652.28 1367.98 1305.27 1519.44 1484.03 1158.84 1280.22 1476.56 1553.9 1418.88 1137.66 1216.41 1347.71 1434.48 41.6177 
6925.23 6250.67 558.706 2408.81 1376.79 291.81 564.634 1325.06 2262.29 2418.84 1722.22 1356.36 2367.68 3805.99 3263.53 1778.92 1309.52 1126.34 1003.31 1260.29 1912.01 2582.63 2228.36 1754.78 1864.88 2063.31 1986.56 2183.72 2161.35 1869.34 1744.42 4.60952 
10313.2 8336.31 219.153 3310.01 1693.2 1017.82 1494.63 1831.67 1536.92 1501.38 2106.2 2032.72 1262.07 1153.69 1661.67 1979.39 2454.6 3552.09 3312.64 1904.19 1799.12 1909.7 1732.43 1447.52 1141.17 1743.32 1755.73 2134.44 2723.13 2133.59 1779.6 97.1057 
28457.2 26407.9 5557.48 8798.83 2777.4 2746.74 2921.43 3974.09 3444.1 3046.34 3716.96 5527.55 6877.13 7677.54 7561.61 4549.97 4470.95 4255.82 4176.11 5254.37 6263.64 6840.92 5019.7 4110.2 2989.92 4138.45 5893.81 3838.35 4991.99 6536.66 5490.01 22.285 
50929.6 52028.3 17994 17751.6 5136.64 9422.76 6748.87 10015.5 14311.3 10087.1 9611.86 11103.9 17280.1 15134.9 11505.6 16095.5 13435.9 8921.56 9902.72 12294.1 11667.6 13102.6 13716.6 9793.89 10582.1 10084 8816.61 13055.1 12757.4 11476.5 10465 56.0126 
76270.3 81537.5 22149.7 20032.3 17047.8 25198.6 21618.9 12459.1 9421.22 11432.3 13330.4 15269.6 16780.3 16602 16237.1 17537 17641.7 18114.9 14257.1 14217 19355.2 22166.6 16678.8 12631.4 14048.8 12253.8 13667 12628.5 10225.4 9155.15 9509.76 227.688 
115357 111207 17879.9 26175.5 28027.4 29323.9 23777.2 39024.6 32334 30885.9 43008.1 25075 24370.2 36515.8 33457.9 26191.2 21404.6 25241.4 21909.4 24117.3 27352.8 25970.5 22451.8 17385.1 16554.9 18405.3 20461.9 24930.6 24907.3 23661.6 24573.6 954.344 
112500 97707.8 20285.4 26311.1 12835.8 13545.5 20452.5 43244.2 22081.5 20481.9 43404.9 27810.3 16223.1 29486.1 20447.7 12543.1 17666.8 26358.6 20637.5 15473.7 22142.8 21647.1 19331.7 21416.9 15976.4 15985.4 18475.6 17892.6 16098.2 19153.7 23801.2 1276.76 
104328 97932.6 13839.7 16408.3 12910.4 15990.7 12264.3 12754.1 10344.9 12517.4 16618.3 32667.7 35642.9 26190.6 19617.7 27169.5 23708.5 13293.7 16237.1 17282.4 14665.9 18232.6 19117.4 22385.9 23062.3 14497.8 14082.4 13716.6 14077.5 14544.9 11828.7 459.283 
150451 134316 11315.5 31579.2 20560.6 23980.5 45076.9 32387.1 18191.9 25467.4 33658.5 35044 35171.5 38993.9 40395.2 50511.7 53111.3 41472.1 28711.5 31435.1 27993 37959.7 50478 45362.3 33935.8 41207.9 40739.4 49512.8 48604.2 35903.3 30197.8 124.755 
172662 154765 17657.3 28014.9 21613.5 21775.7 33100.8 51828.6 48906.9 38970.3 36779.2 26117 19249.7 18693.5 47662.3 66902.2 26841.9 14892.1 30103.8 36930.6 19141.5 20345.8 38046.8 37887.4 20417.6 25519.1 23764.1 25407.6 19245.5 28353.9 25150.5 22.1303 
97870.8 70121.5 18496 18882 18658.9 36276.3 22682.4 14523.2 22291.8 21751.4 18708.2 10289.8 17555.9 15770.5 10274 12777.6 15642.3 24406 20930.1 13506.8 26811.6 30722.4 18508.8 25166.4 23398 31523.3 16310.6 14142.4 17905.8 23228.4 19110 206.442 
106895 76193.9 32048.4 31914.8 24154.4 58788.7 40586.1 25272.9 20358.6 39078.5 25586.7 29373.5 28543.8 23073.9 28993.9 30509.2 27650.3 42116.7 41475.1 34482.4 28423.1 26305.1 37905.4 46409.3 29842.9 22602.6 31820.9 44647.3 30919 31150.5 40283.4 1168.05 
110054 117029 59223.7 45508.3 34261.5 37161.3 17824.2 22048.1 27322.3 28776 21080.2 20734.5 27030.1 25667.2 18042.2 32198.8 26688.1 24623.3 31426.4 32478.1 21286.2 16348.5 24431.4 25795.3 20545.7 32672.6 35039.5 29841.4 28351.7 30023.9 25184.2 478.571 
155972 225717 127469 78579.1 43055 32368.8 41324.1 43503.5 42339.5 39289.7 26537.9 25414.6 21008 21862.9 23082.1 27518.5 34992 26760.7 20728.6 32118.5 40097.8 43008.1 49123.9 46466.5 32704.2 31723.5 36574.6 36473.8 29732.1 35067.9 39399.5 3067.23 
156864 236025 133409 91513.4 57538.8 46021.4 49735.7 27997.9 27336.8 50995.7 48346.2 35143.2 30042.2 28570.9 26422.7 22353.3 25702.2 34108.2 35059.5 37702.2 44111.4 50991.7 65494.5 82688.4 88328.6 85737.1 79783.2 70612.3 62840.6 68213.4 71880.6 4687.58 
152427 250000 147989 98998.5 83579.7 62064.6 48330.7 28818 29697.5 53410.7 58406.2 49217.5 47462.5 47659.9 38607.8 31770.2 39402.2 45704.9 39796 37310.1 42531.1 46941 51969.1 62508.7 72001.3 70994.6 64425.4 58271.9 49121.3 46642.1 56735.9 5019.48 
114019 196918 121818 72337 58986.6 41050.3 23889.5 9632.43 7121.27 17776.3 23692.4 25003.9 29463.6 28454.7 21849.3 19351.8 19543.3 19503.3 21508.7 24321 21289 18675.5 22415.4 28023.5 32358 37221.5 37658.2 27980.1 24130.1 35428.2 46181.4 3307.51 
99110.2 186192 132665 83791.5 73607.4 64189.8 51193.5 38362 30950.7 31281.6 31258.3 28794.3 28872.1 30964.7 31862.6 33249.1 38032.4 45869.4 50843.6 49193.9 48317.3 51068.3 51993.4 52895.8 54505.9 53322.1 52179.4 50220.4 41814.6 31800.1 32217.3 6053.79 
87901.4 174426 132633 79278.7 65702.2 57041.6 42773.6 34072.3 30899.8 31910 35100 38120.8 41637.3 42223 37822.6 33958.9 33684.9 37169.8 45948.1 54595.6 53778.5 44772.2 35918.9 30746.6 28636.5 28465.8 29471.3 32774.4 37957.1 41318.8 43533.8 6074.29 
81821.2 175358 154223 101039 75681.9 60800 41588.6 32684.6 35666.1 36535 33037.3 31797 33105 32814.2 29494 25214.5 22974.3 23934.2 26938.8 30580 34777.9 39866.7 44173.6 44698.1 41882.5 39915.9 40562.5 42054.4 44083.6 47184 49201.7 2459.96 
27233.4 57115.4 45557.8 23104.9 16306.7 16740.9 10362.6 5056.77 11974.9 22737.4 23498.9 16975.3 14901.1 18289.2 18613 13270.6 8838.94 10289.1 14607.1 16228 14445.1 13124.6 14621.1 16635 16016.7 13708.3 13398.7 16001.5 18609.3 19472.4 20794 1229.9 
25673.6 54696.7 45562.5 25618 20235 23106.2 20172 15629.3 18487.1 23255.2 20978 16001.3 16958.4 21299.5 21650.6 19233.8 19553.9 21693.7 22034 22558.1 26237.3 29456.4 27463.6 22713.8 20215.3 19086.2 16039.7 13019.9 13965.9 17345.1 18754.5 3931.2 
10894.2 25196.9 25019.9 16956 10712.7 8578.45 7390.47 5807.67 5141.1 5461.79 5259.74 4338.9 3924.32 4225.4 4391.28 4500.42 5222.93 5853.94 5133.38 3863.27 4111.57 5615.8 5710.88 3584.61 1907.02 2810.04 4649.67 4820.69 3921.79 4509.23 6705.82 1944.74 
4854.94 11182.4 11146.1 8182.51 6868.9 7437.05 7136.77 5300.86 4362.28 5905.52 8089.64 8104.35 5779.18 3384.21 2578.15 2961.51 3451.98 3715.06 3852.19 3701.37 3160.43 2656.14 2610.58 2730.88 2435.33 1808.42 1505.42 1728.14 1932.28 1732.98 1524.49 200.878 
998.22 2746.83 3902.95 4426.18 4466.13 4165.02 3608.19 2914.79 2276.89 1858.82 1685.33 1668.62 1733.76 1876.48 2092.72 2298.63 2365.21 2232.03 1961.54 1676.44 1461.69 1331.48 1262.59 1226.83 1190.41 1111.91 966.633 773.806 588.635 462.999 415.835 17.5075 
779.939 2260.85 3516.39 4432.72 4943.95 5039.32 4761.3 4195.41 3454.01 2657.29 1915.02 1312.04 899.493 692.634 674.594 804.583 1028.4 1289.05 1535.68 1729.68 1847.53 1880.79 1833.81 1720.27 1559.31 1371.94 1178.27 995.373 836.114 708.492 615.591 103.458 
911.082 2670.78 4249.61 5547.96 6494.78 7053.76 7225.22 7043.56 6570.91 5888.06 5083.94 4245.17 3447.05 2747.03 2181.14 1763.62 1489.27 1337.89 1279.9 1282.06 1312.67 1345.4 1361.68 1351.36 1312.09 1247.73 1166.17 1077.1 990.141 913.364 852.459 121.66 
];

change60 = [0.0795775 0.238732 0.397887 0.557042 0.716197 0.875352 1.03451 1.19366 1.35282 1.51197 1.67113 1.83028 1.98944 2.14859 2.30775 2.4669 2.62606 2.78521 2.94437 3.10352 3.26268 3.42183 3.58099 3.74014 3.8993 4.05845 4.21761 4.37676 4.53592 4.69507 4.85423 6.09959e-027 
9.62887 15.0706 9.85578 12.5895 13.4945 8.39229 7.84119 10.0993 11.9729 13.7054 11.7316 8.78113 10.7104 17.3568 26.3859 34.7699 38.3114 33.8398 27.6912 31.2331 36.1944 31.1421 30.4959 39.3304 40.9106 35.366 32.0396 34.6943 46.0984 53.8533 56.8669 0.00901218 
45.8366 64.6391 46.9024 46.7512 33.1498 40.2239 43.109 35.924 28.0586 25.8084 33.4802 26.7922 27.9163 49.3825 72.8511 68.5895 54.6743 52.0823 50.6465 81.7211 90.8522 63.673 69.016 84.423 80.2974 57.7767 76.3802 121.082 121.43 124.404 128.639 3.26415 
548.209 656.489 264.4 243.76 318.374 485.97 425.732 310.826 229.845 236.503 303.313 313.427 358.218 375.898 310.875 255.75 301.997 363.492 409.061 422.662 350.593 320.189 300.311 269.43 261.275 319.359 347.671 420.268 412.659 387.63 495.769 12.6268 
4136.76 4909.87 2012.02 918.144 844.255 1688.22 1623.63 1637.1 1035.1 1303.34 2073.7 1988.82 1358.41 1092.58 1218.15 1546.07 1654.27 1684.32 1302.17 1193.77 1385.49 1509.95 1199.85 1305.64 1567.71 1422.65 1411.49 1208.99 1038.88 1230.84 1544.24 56.7827 
6925.23 5063.85 620.091 1635 1263.46 1114.05 677.931 1961.45 2609.3 2362.72 1559.66 1500.22 2654.37 3380.43 2815 1949.79 1038.15 1854.49 1166.87 1395.48 1976.72 2546.38 2291.69 2303.13 1896.33 1606.35 1903.97 2440.02 1961.83 1695.93 1863.61 49.4128 
10313.2 7548.44 487.796 3340.12 2080.87 1650.52 1147.17 1864.59 1070.86 1566.33 2650.89 2120.09 1381.53 1232.36 1434.89 1883.62 2609.19 3347.84 2943.27 1948.83 1382.85 2202.47 1656.97 1498.36 1546.34 1506.39 1715.76 2291.33 2485.61 2173.73 1853.3 110.396 
28457.2 19057.6 5512.04 7924.31 3410.17 6038.92 3448.15 5305.39 2539.74 3090.77 3594.22 5486.99 5037.24 8779.83 6255.7 4710.54 5736.22 4500.62 4885.79 4903.15 6430.18 7363.28 4735.75 4488.88 4076.83 4831.1 6203.82 5071.71 5321.72 5381.29 4780.98 9.68271 
50929.6 30789.7 18013.6 14160.7 4617.18 13482.5 5182.16 9701.21 12171 11836.3 9024.48 16589.2 14644 9258.78 11294.5 11139.9 11085.2 13098.4 8210.65 14336.8 13107.3 12412.6 12899.5 11143.2 9925.78 10268.1 9599.42 10727.3 13333.2 10698.2 10307.5 184.897 
76270.3 54419.8 21988.8 15217.8 15055.6 35906.5 24774.5 17225.7 10821.7 21325.8 11700.6 16373.4 18238.9 13136.5 15779.8 20044.3 18198.6 17054.1 14309.3 16335.6 17182.2 22280.9 14545.5 16141.8 14899.6 11739.4 13383.8 12493.5 11404.3 10770 7626.54 234.779 
115357 113155 15690.5 14034.1 30473.1 50667.2 23432.6 21100 37197.8 43519 37689.5 31590.9 32102.9 27316.6 27560.7 29859.8 24245.3 19671.6 20422.4 28567 24855 21979 22704.3 17964.4 15365.5 15878.3 21886.8 25490.5 24812.8 24263.5 24434.9 752.503 
112500 118392 19554.9 5506.64 13249 37512.2 18597.1 14034.8 28445.5 41198 40399.8 26447.3 20799.6 16390.5 18908.2 18644.3 15923.7 19835.1 20559.1 21591.6 19407.5 19419.3 18660.1 16481.1 16769 18313.9 19613.3 15607.5 16078.9 21419.4 22912.8 578.837 
104328 97181.5 13876.6 21825.1 14233.4 21334.5 10412 8195.99 12082.5 17280.1 16982.8 25178.2 41181.3 28045.6 21206.6 27394.2 16879.6 12711.6 17028.6 15145.6 15549 18401.4 19005.4 21042.3 21068.4 15599.8 14670.8 16935.4 13865.2 10670.7 13026.6 988.402 
150451 132906 14551.9 34208.9 15505.6 38576.4 41160.6 21741.5 22000.7 30561.1 35613.5 43421.7 40024.9 35401.6 38804.5 38051.8 50066.5 46286.5 29268.6 30617.5 30420.4 34817.6 44282.2 43855.8 48918.3 43242.8 33900.7 45408.1 50227.3 40611.4 30789.2 304.754 
172662 168552 23151.9 44479.3 23712.8 17874 37248.8 44825.3 46646.2 36997.6 32650 39853.8 24424.6 24474.8 51430 49408.2 27528.5 24771.1 33748.1 31991.8 18704.2 24498.6 39462.7 29767.6 18071.1 26425 26223.7 21739.9 23705.6 27206 25627.4 1487.55 
97870.8 87203.8 24727 27201.4 24579.5 38496.1 19708.1 22047.7 26604.6 15935.6 16061.8 12679.7 17947.7 19344.1 12032.3 19042.2 21199.8 28723.3 19946.4 16459 32191.5 32495.2 25160.7 23747.7 31744.2 34714 17415 14018.6 19456.2 18535.3 12659.4 539.386 
106895 91454.7 42481.2 56583.9 38008.9 67799.8 39655.6 23904.9 34357.5 38211.3 29233 31863.5 38615.2 30420.9 40169.4 37914.2 31418.8 50232.5 45603.6 37043.4 24124.1 31536.8 46899.9 34359.5 27359.5 21945.2 27495.9 30901.6 21015.2 23409 30831.4 206.15 
110054 104061 61023.3 75190 41145.8 37251.1 18559 27117.2 30664.2 25665.3 26077.6 25801 29594.8 29655 25728.9 45098.4 30870.1 21277.2 35915.8 29011.1 19672.2 18821.3 23979.1 25001.2 26741.4 33349.9 35322.7 31208.5 31534.8 35453.3 28715.3 124.422 
155972 212755 127227 86918.2 44159.3 41144.6 38767.2 41099.6 43479.8 36159 28422.8 30540.6 21899.8 19379.8 24389.5 32743.5 31940 24160.4 23561.4 27314 36094.7 49330.6 52053.5 40482.5 29987.3 36482.6 37036.4 32577.6 30978.9 27102.7 39370.2 6064.3 
156864 213623 133506 78670.1 57157.8 88922.9 49228.4 30159 27496.2 34439.6 48414.9 44004 29672.2 22934.1 26757.8 27049.3 27792.3 31952 32146.7 35256.7 43938.3 55322.1 68300.9 82348 85638.6 83151 79288.9 70726.6 64588.4 62995.8 71795.4 7153.37 
152427 214581 148089 81388.6 82309.4 115335 48537.4 39055.3 32496.3 37207.5 54320.1 55332.2 49900 41868.9 39441.3 36398.9 38228.3 38102.7 40837.3 42264.1 41035.2 47125.9 52615.6 64643.7 69191.5 69549.3 66753.8 55966.9 49823.4 48892.3 57151.7 4338.99 
114019 163124 121193 66761.5 56970.2 74730.9 23976.1 16088.6 8438.91 12703.9 23151.6 27932.4 31258.8 25783.4 21301.8 19799.3 18427.1 17872.4 18341.7 23408.2 23251.1 21323.4 22716.6 27190.1 34248.9 38922.1 37025.6 27055.6 23818.1 27436.1 46180.3 5953.96 
99110.2 152595 132274 84086.4 72191.2 98938.6 49834.7 37277.7 33376.2 28310.9 30739.2 30740 30151.5 30438.6 30992.5 30481.9 37130.4 50362.8 54181.1 48212.5 47737.7 47069.7 48588.9 54782.1 50874.7 52258.7 54729.2 51582.7 44105.1 35187.7 35193.2 7616.24 
87901.4 155873 133206 79955.7 66734.5 68742.2 41804.3 34868.3 30552.5 34587 35955.4 40667.4 38395.7 39652.4 38792 36249.1 35775.6 38446 46474.5 50901.9 53625.3 43571.2 35925.3 32898.5 27994.4 27928.5 29017.7 34243.3 37838.5 40407.4 42662.5 6564.55 
81821.2 166377 154340 101190 75606.8 64210 41380.1 32365.6 36138.2 38732.4 34298.5 33583.3 32155.4 34473.9 28461.5 25021.4 22866.6 22904.8 26211.4 30204.9 34687.7 41797.5 45838.1 46532.4 42367.8 39654.4 39455.3 40839.2 43664.3 46760.6 49871.8 2230.88 
27233.4 56217.6 45557.4 22907.5 16272 16269.2 10299.1 4934.44 11859.4 22795.4 23388.8 17586.7 14875.7 18983.2 18642.8 13427.2 8981.08 10322.7 14761.7 16148.6 14533.4 13200.2 14600.2 16726.9 15934.1 13544.8 13302.3 15866.4 18568.9 19265.9 20828.1 1140.19 
25673.6 54696.7 45562.5 25618 20235 23106.2 20172 15629.3 18487.1 23255.2 20978 16001.3 16958.4 21299.5 21650.6 19233.8 19553.9 21693.7 22034 22558.1 26237.3 29456.4 27463.6 22713.8 20215.3 19086.2 16039.7 13019.9 13965.9 17345.1 18754.5 3931.2 
10894.2 25196.9 25019.9 16956 10712.7 8578.45 7390.47 5807.66 5141.1 5461.79 5259.74 4338.9 3924.32 4225.4 4391.28 4500.42 5222.93 5853.95 5133.38 3863.27 4111.57 5615.8 5710.88 3584.61 1907.02 2810.04 4649.67 4820.7 3921.79 4509.23 6705.82 1944.74 
4854.94 11182.4 11146.1 8182.51 6868.9 7437.05 7136.77 5300.85 4362.28 5905.52 8089.64 8104.35 5779.18 3384.21 2578.15 2961.51 3451.98 3715.06 3852.19 3701.37 3160.43 2656.14 2610.58 2730.88 2435.32 1808.42 1505.42 1728.14 1932.28 1732.98 1524.49 200.879 
998.22 2746.83 3902.95 4426.18 4466.13 4165.02 3608.19 2914.79 2276.89 1858.82 1685.33 1668.62 1733.76 1876.48 2092.72 2298.63 2365.21 2232.03 1961.54 1676.44 1461.7 1331.48 1262.59 1226.83 1190.41 1111.91 966.633 773.806 588.635 462.999 415.835 17.5075 
779.939 2260.85 3516.39 4432.72 4943.95 5039.32 4761.3 4195.41 3454.01 2657.29 1915.02 1312.04 899.493 692.635 674.594 804.583 1028.4 1289.05 1535.68 1729.68 1847.53 1880.79 1833.81 1720.27 1559.31 1371.94 1178.27 995.373 836.115 708.493 615.592 103.458 
911.082 2670.78 4249.61 5547.96 6494.78 7053.76 7225.22 7043.56 6570.91 5888.06 5083.94 4245.17 3447.05 2747.03 2181.14 1763.62 1489.27 1337.89 1279.9 1282.06 1312.67 1345.4 1361.68 1351.36 1312.09 1247.73 1166.17 1077.1 990.141 913.364 852.458 121.66 
];

figure,
bar3(origin);view([60,60,60]);set(gcf,'renderer','zbuffer');set(gca, 'ZLim', [0.0 300000.0]);
figure,
bar3(change30);view([60,60,60]);set(gcf,'renderer','zbuffer');set(gca, 'ZLim', [0.0 300000.0]);
figure,
bar3(change60);view([60,60,60]);set(gcf,'renderer','zbuffer');set(gca, 'ZLim', [0.0 300000.0]);


% 
% figure,
% bar3(theta_change);view([60,30,60]);set(gcf,'renderer','zbuffer');%set(gca, 'ZLim', [0.0 15.0]);
% figure,
% bar3(positive_N);view([60,-60,60]);set(gcf,'renderer','zbuffer');%set(gca, 'ZLim', [0.0 15.0]);
% figure,
% bar3(theta45);view([60,-60,60]);set(gcf,'renderer','zbuffer');%set(gca, 'ZLim', [0.0 15.0]);
% figure,
% bar3(theta45_negative);view([60,-60,60]);set(gcf,'renderer','zbuffer');%set(gca, 'ZLim', [0.0 15.0]);
% figure,
% bar3(phi_change);view([60,-60,60]);set(gcf,'renderer','zbuffer');%set(gca, 'ZLim', [0.0 15.0]);
% figure,
% bar3(image_negative_N_negative);view([60,-60,60]);set(gcf,'renderer','zbuffer');%set(gca, 'ZLim', [0.0 15.0]);
% figure,
% bar3(image_negative_N_positive);view([60,-60,60]);set(gcf,'renderer','zbuffer');%set(gca, 'ZLim', [0.0 15.0]);
% 
% 
