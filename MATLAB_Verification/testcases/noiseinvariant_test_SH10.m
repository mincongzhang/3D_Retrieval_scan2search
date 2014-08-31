clc
clear
close all

input = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
1.97466 1.33504 2.19498 2.5614 0.960229 1.9545 1.53613 2.58094 2.40028 3.03546 4.51238 3.68061 4.36602 4.2566 3.06898 4.62276 4.41037 3.96596 4.318 3.98392 4.52632 4.95044 5.12135 5.52039 5.5891 5.41027 5.90766 5.74076 5.06857 5.51435 5.08592 0.432864 
6.20609 5.03928 3.53252 5.59791 3.28215 2.55245 3.03794 3.54805 2.64104 3.29799 4.42618 5.89042 4.96935 7.32409 9.34363 8.27167 7.6889 6.07201 7.51045 7.58054 8.69759 9.59679 8.12537 9.90873 11.0168 11.6812 9.94862 8.50021 9.12927 9.07464 9.10141 1.30235 
9.59122 8.80136 1.7391 4.71679 7.11193 7.18364 4.92923 4.27678 3.11404 4.80629 5.16528 5.86009 5.33892 5.83337 7.80321 8.53173 8.42159 10.4394 8.77859 9.12327 11.276 10.3767 11.794 9.68281 9.53431 10.6112 11.3704 13.0876 16.0708 17.8918 15.7787 0.890673 
15.7973 16.1863 7.3961 4.26469 6.92845 9.48142 10.2436 9.19929 8.81771 8.03213 6.24727 7.5589 4.97972 5.40829 5.89877 9.37319 8.60505 9.19191 8.97022 10.6313 14.8394 16.0648 13.892 16.3366 14.207 13.4976 17.7948 16.0095 18.6587 17.3797 14.2217 1.86708 
22.0034 24.5753 16.066 12.4571 6.97762 10.3263 15.6817 18.514 19.0364 17.5731 14.723 12.1376 10.0269 10.4013 11.643 11.5247 11.2335 12.1822 10.9038 10.5062 11.1347 13.0421 13.0178 12.1147 13.4218 15.2243 17.3484 18.3822 17.0236 17.4608 17.658 0.906264 
19.7466 19.8179 17.6799 18.6289 10.8797 10.9567 11.3762 13.9552 17.2259 17.5834 18.5263 18.8748 16.8434 16.3404 14.2818 13.0575 11.7964 12.2345 13.7599 14.0536 14.0616 13.8698 10.9273 11.8403 13.7923 15.3058 16.1753 15.6813 15.6454 15.9312 14.5645 2.10139 
11.0017 10.3794 11.5438 10.8874 9.76024 10.4853 9.25699 10.952 9.3992 12.054 12.4226 10.6597 12.06 12.2574 11.61 11.1869 10.6843 10.6069 9.86771 10.6184 10.1238 10.9916 11.3111 10.6361 10.3538 11.7312 10.064 11.6434 12.1221 12.7292 12.2828 2.38438 
11.5659 14.3182 14.3685 14.862 13.8662 14.1255 11.0874 7.58653 6.43179 8.47793 9.08453 8.60502 11.1917 10.8202 10.6614 12.2533 12.044 12.1952 12.7461 12.5222 12.2349 11.7673 11.1928 10.6905 10.2526 9.60284 9.5894 10.725 10.9339 10.2421 9.26005 2.99801 
8.74494 8.90119 10.0631 11.4903 8.94458 11.3731 10.7989 9.12785 8.5498 8.97061 10.5368 9.68667 9.18983 9.34976 9.29663 9.68136 10.7737 11.2248 12.1859 12.4409 11.1135 10.808 9.97572 9.52111 8.72007 8.50972 9.92385 10.1845 11.8455 12.2174 11.4328 2.06671 
11.2838 10.8823 13.1769 13.2148 11.5802 16.428 12.8116 10.88 8.96872 8.92636 9.02821 7.97855 10.0037 8.72571 9.81095 10.7639 12.6509 13.4271 15.3529 15.2498 13.2726 12.9881 12.4435 14.0985 12.5582 10.8669 10.2388 10.7208 11.8565 10.7909 10.8448 1.20679 
442.889 509.885 477.509 451.561 403.084 362.509 289.769 253.142 181.312 174.471 111.178 110.89 83.6724 90.1666 83.2013 84.4548 99.9307 85.7744 93.4926 90.536 94.9078 109.429 96.7599 133.756 105.273 108.629 102.874 103.999 111.297 118.275 109.34 16.3265 
737.396 776.968 586.254 502.454 220.44 257.903 325.341 419.797 484.748 469.895 401.284 330.268 201.831 186.189 159.431 208.731 226.196 217.453 184.651 156.245 110.356 132.075 144.727 151.946 173.099 187.262 189.616 179.616 143.941 115.977 124.714 10.8935 
190.132 188.942 70.9729 75.3461 136.674 181.657 194.306 165.537 93.6651 64.6856 92.7937 124.571 131.897 111.463 78.9902 71.0214 70.18 81.4231 82.3347 75.0107 63.6314 62.0188 58.8626 62.3877 65.2226 65.367 67.4151 69.9585 56.9038 60.0297 66.1073 16.2058 
124.968 108.586 15.3467 79.2295 139.883 131.859 77.9897 53.659 77.1345 106.32 98.5553 72.1409 30.1547 63.5181 70.9443 64.3708 40.8014 42.4974 37.8162 45.6765 44.99 35.1952 32.5355 36.6713 38.846 36.1546 37.2463 39.5962 36.2159 35.9344 40.2335 4.07343 
173.206 142.976 102.936 190.835 206.01 124.249 81.4788 177.356 220.564 133.714 103.715 183.993 193.687 135.537 109.019 171.103 188.231 136.568 99.9131 162.191 173.549 130.605 103.062 153.063 165.916 123.228 93.5823 141.891 152.729 119.793 94.7398 25.0584 
179.13 135.57 155.161 232.463 168.415 70.9053 158.137 220.562 168.221 73.9187 207.025 233.184 114.65 129.261 219.668 187.172 89.6726 142.396 211.622 171.759 112.038 169.177 195.168 139.906 90.6582 166.347 174.316 113.588 102.654 165.481 148.857 11.3047 
152.331 100.035 173.509 208.609 66.1303 92.0864 172.968 122.208 54.5394 123.219 140.285 87.0768 89.9989 146.288 112.994 63.5917 137.528 131.259 51.7215 111.303 110.982 70.8058 79.2657 89.1334 61.2047 66.6748 66.788 69.9385 61.3547 58.6414 56.7548 11.0929 
205.365 118.533 257.977 281.586 56.3756 192.477 290.161 154.804 183.089 276.751 129.779 150.431 267.537 181.742 108.857 232.478 174.395 94.4672 209.91 189.841 64.9423 173.097 187.805 84.9955 149.82 189.567 87.5915 127.466 183.209 111.402 113.175 43.4092 
110.299 64.4487 145.553 157.091 17.6856 113.758 153.36 60.6883 124.361 153.553 46.3113 105.786 147.741 73.4252 98.3465 136.129 68.3047 88.5564 134.825 93.7443 80.4009 121.538 74.1355 61.9806 110.177 82.2185 64.3771 101.423 72.1499 55.7771 87.2879 10.8736 
156.845 84.915 216.705 215.361 17.4064 177.622 208.706 79.022 204.841 207.732 33.8076 176.569 199.222 83.4435 180.644 191.782 41.4143 155.776 182.126 81.9747 150.622 172.487 52.8039 127.962 165.043 79.4824 122.291 157.024 64.8386 108.03 149.122 4.77337 
58.1115 27.9442 92.3936 78.1595 36.487 85.1553 48.2714 41.0901 90.3222 47.8713 65.5501 83.1838 30.4537 71.7321 68.6838 51.9071 69.5338 68.1976 41.1055 68.516 52.3558 48.5589 62.0139 44.6531 41.6218 45.9703 29.0741 36.8092 41.4068 36.9474 39.099 18.773 
41.75 17.7327 73.2148 50.7013 47.1837 65.862 11.9621 50.0954 57.1298 20.7303 71.8491 43.8993 50.3089 64.4236 29.9015 59.3389 53.3184 46.7421 61.3334 52.5996 41.572 57.9323 23.5023 46.1274 41.5411 31.1811 48.2656 38.1273 36.0202 46.1411 27.7108 9.73 
45.4173 16.3042 83.4572 51.6751 64.0782 74.1762 12.6963 68.6555 45.7056 37.2267 79.4869 25.3568 77.417 59.6203 49.6038 76.3398 42.9615 70.3151 66.229 57.9157 70.5185 58.5602 48.2646 62.0728 23.3731 53.7556 41.7052 40.182 56.0563 41.2887 48.8084 14.0776 
67.7028 21.3843 130.888 66.3958 119.393 101.834 61.4939 108.489 24.8093 81.4826 87.231 35.3527 120.611 46.3881 117.399 89.452 96.4272 110.175 92.5693 106.047 107.15 90.3278 108.966 81.6862 85.6398 82.6652 49.9822 80.0992 49.1754 69.7177 78.9364 18.0251 
124.686 34.8317 248.766 117.227 249.549 190.454 170.876 223.012 52.3194 201.378 74.7001 134.005 155.314 60.6874 179.326 88.5677 151.068 143.057 105.893 162.028 105.76 148.033 138.106 126.623 148.402 125.529 125.067 137.4 88.604 138.562 84.6518 44.5197 
149.792 41.9021 301.387 138.492 309.062 229.199 222.015 277.818 85.2738 265.452 97.5094 195.724 211.261 104.363 263.432 109.316 240.757 188.818 163.083 235.289 105.135 229.758 161.501 187.714 227.589 155.719 241.405 175.762 195.732 212.668 114.488 88.5148 
95.9122 29.5409 192.309 99.9564 195.056 166.056 135.919 201.388 50.7254 191.22 87.968 135.976 167.346 55.7098 205.704 68.4957 192.837 144.845 140.012 189.629 89.4532 189.998 112.357 151.528 161.743 100.288 181.741 92.8154 161.167 131.576 112.084 76.12 
56.9831 13.1266 114.639 44.9398 117.852 76.0868 87.5662 96.7224 54.2699 104.056 75.1172 104.571 115.249 110.457 134.175 125.408 122.52 139.268 87.0445 140.849 59.8648 127.281 86.543 107.086 126.578 98.3106 147.326 108.71 142.226 122.574 119.022 4.41722 
3.38514 1.46908 6.92557 4.09238 7.45283 6.84079 6.09772 8.78977 4.16986 9.37113 4.58393 8.50296 7.24635 6.80613 9.64881 5.97124 10.8231 7.4244 10.6252 9.80087 9.51595 11.5174 8.50804 11.8978 8.63934 10.9078 9.79375 9.08979 10.9498 7.64227 11.3376 0.759219 
0.282095 0.488603 0.630783 0.746353 0.846284 0.935603 1.01711 1.09255 1.16311 1.22962 1.29272 1.35288 1.41047 1.46581 1.51913 1.57064 1.62051 1.6689 1.71592 1.76168 1.80629 1.84982 1.89235 1.93394 1.97466 2.01456 2.05368 2.09207 2.12977 2.16681 2.20323 0.7107 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
];


database_data = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
2.82095 2.71888 1.05078 2.34252 2.4494 1.93877 1.30106 2.91171 4.13131 3.67493 3.86216 4.48139 4.64754 5.11155 5.59118 5.91168 5.30614 4.5713 4.62461 4.01131 4.72296 6.09915 6.68588 6.99043 6.22811 6.40808 7.1605 5.90295 6.28485 7.47224 6.88133 0.278405 
4.51352 4.37583 3.03351 2.50936 1.59078 3.00395 4.3151 5.05278 5.26247 5.26657 3.76724 4.80564 3.43251 5.11793 6.87887 6.76232 7.27575 5.79817 6.36587 6.19281 7.1108 7.72719 7.45951 8.4136 7.60302 8.00236 8.91435 7.90093 5.99562 7.82002 8.3511 0.00862025 
5.92399 5.81178 3.60897 4.04714 4.20161 2.64288 3.62816 4.4592 3.21119 4.83166 7.67804 8.29567 5.85335 5.05609 6.74409 7.47744 4.40506 5.82373 5.51415 7.30521 8.44544 9.69722 10.5707 7.9374 9.14389 8.81561 9.65227 8.97797 8.37823 9.09039 11.1833 1.10513 
15.2331 11.8794 4.409 6.86052 7.89682 12.5958 8.57459 4.95727 8.24852 8.87362 4.96921 5.38251 5.65659 8.51927 9.96464 11.1853 11.3754 9.9868 6.60133 6.39966 5.81032 9.29741 7.41282 11.1453 9.03586 17.1439 20.8918 17.311 21.1726 17.0305 19.3898 8.33987 
19.1824 17.1675 4.1184 7.58657 4.36668 6.14147 11.0814 14.2718 12.2641 7.45566 8.95842 13.9725 15.1209 11.7465 5.30262 7.195 5.89326 8.44356 15.4626 17.5792 15.984 13.0326 7.97398 7.54947 7.33678 9.52504 9.48532 10.7815 9.10984 15.9249 13.277 2.35743 
24.2602 23.7696 12.5871 8.16758 5.15648 10.4993 4.79862 8.23136 12.7306 13.961 11.8864 8.82278 9.77237 14.9748 19.8805 18.4565 14.1027 11.3455 7.48193 9.31975 8.92492 11.9896 16.4783 19.4107 21.3532 19.789 13.7093 11.9607 10.0793 10.163 10.2249 1.46704 
54.4443 56.0591 38.4171 30.6995 6.97713 14.2766 17.6827 20.6405 15.4517 14.4034 9.31614 9.79262 10.4929 13.0067 13.6952 17.5649 25.2727 29.0543 28.2986 25.2539 21.4558 17.6626 13.9089 14.6671 14.8021 18.4755 22.7607 23.6291 23.3661 22.7053 21.4427 3.70897 
17.772 20.5784 13.292 4.922 9.97846 16.6693 21.0852 20.6221 18.6683 13.8385 9.7517 11.3038 15.0259 17.4383 17.5844 14.6572 9.37575 9.82657 11.4312 14.0781 15.4751 15.6311 13.0548 10.2531 9.40928 11.4021 14.5444 14.8798 12.5471 12.2542 11.1781 2.08625 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
16.9257 23.3661 22.9981 22.9994 20.6229 18.957 19.212 17.3767 16.7899 15.9802 12.5982 11.9173 11.7453 9.87955 10.8791 16.3587 21.6005 20.174 14.2468 13.819 11.9448 8.9464 15.8721 18.7762 18.7476 19.9404 16.6206 13.5682 17.0158 16.3037 12.8974 0.935712 
477.586 523.202 475.442 438.906 336.285 291.417 178.99 158.023 69.5971 89.9077 58.3927 81.3111 58.7578 61.3274 58.6275 71.6584 61.8802 67.1865 63.7385 61.0102 42.2969 45.7853 37.3472 40.5747 45.1858 49.7869 46.7953 55.6706 56.7031 72.3253 75.4857 17.8441 
746.987 783.044 603.943 503.344 223.924 211.449 179.893 266.113 286.482 278.012 217.106 187.199 97.0347 106.349 76.7011 97.5246 94.6978 134.054 100.831 96.6736 79.2212 77.6861 63.7828 74.2612 70.4888 77.9998 76.8973 83.1164 81.7452 88.2707 95.8849 30.8965 
223.983 229.359 92.0973 82.9437 146.996 186.389 189.431 167.015 112.762 80.1555 65.483 80.525 91.7618 95.8177 88.1142 72.6641 49.3139 48.2148 69.545 80.3447 67.6543 52.469 34.9341 47.2559 59.6388 60.4925 48.7208 44.0517 43.299 43.6631 46.6114 13.5536 
136.534 114.965 28.9764 104.426 156.241 133.74 39.7896 58.3374 107.508 108.455 56.5749 39.7119 59.6326 68.3858 40.7271 31.9838 39.5857 47.3909 32.9513 31.5089 43.0352 44.4892 28.2499 38.3248 47.6462 42.9921 27.715 37.9842 42.9408 42.1723 33.616 5.71697 
228.215 195.579 118.05 243.977 282.564 179.806 82.7258 202.747 253.522 171.113 101.773 187.285 208.709 153.746 110.602 154.993 164.74 130.344 113.239 147.467 146.96 121.163 110.728 140.671 143.35 112.908 110.694 144.464 143.196 117.766 106.132 21.6168 
209.314 172.052 119.799 236.202 229.047 114.865 112.273 197.87 166.611 78.0778 121.993 149.89 89.8678 87.4116 132.483 106.971 72.389 129.536 148.086 98.5751 99.3752 164.829 166.096 105.892 101.634 178.04 176.87 108.631 97.3841 168.099 163.392 28.3182 
93.9376 54.8687 109.103 117.176 27.2189 70.4875 97.0724 53.025 51.9329 68.5456 25.7711 35.2019 39.6819 29.1721 21.5288 27.778 23.0244 36.5072 41.4403 33.6934 45.9814 66.6631 38.7948 43.0134 80.2 56.2928 44.0933 72.6746 47.598 46.9859 62.0915 3.58 
256.424 151.002 322.836 357.701 60.5722 239.05 358.171 178.165 229.602 340.183 153.236 176.076 325.726 212.767 122.346 286.198 206.869 109.446 264.172 233.184 57.968 226.085 230.506 89.6544 206.548 239.626 74.7112 182.985 237.084 114.753 172.384 69.8341 
84.0642 44.9738 117.234 109.928 21.1596 84.4053 68.3399 23.4421 50.625 60.3471 52.4181 50.5299 65.3013 72.97 48.3869 60.9036 92.4527 60.174 61.6467 88.8204 58.0224 54.7624 94.9387 71.6979 62.2875 87.2709 82.3518 61.6455 84.4249 87.5477 63.788 16.1473 
142.458 78.0444 205.052 199.268 39.7859 174.231 145.512 61.9985 157.333 135.798 35.6628 117.346 89.1739 51.1454 65.8996 79.3393 46.0046 45.6431 68.7108 81.0142 40.6248 69.902 114.87 69.1601 90.0229 135.022 71.4306 94.364 151.951 87.8786 111.19 55.9282 
105.221 52.6622 162.562 144.012 56.4601 149.037 94.6353 73.6599 144.68 96.6155 75.4465 121.495 47.0425 68.4444 93.4581 61.2055 53.5411 77.9237 40.9974 42.0468 60.7758 53.8328 32.2371 56.174 62.8426 32.3937 67.4356 71.748 37.3237 71.3435 78.5095 3.66224 
32.4409 11.9727 60.4571 37.1422 48.8893 54.4629 14.6261 53.1855 27.2506 33.4761 54.1948 17.3825 57.9181 38.7827 39.5229 54.3766 18.691 52.6496 35.7438 40.3117 51.5163 37.7179 49.1348 48.4801 33.4772 53.1197 28.2101 45.3053 41.9502 33.1885 48.2721 6.77896 
39.7754 13.3788 76.1585 42.6229 67.1679 65.0904 29.6817 68.2257 20.4097 49.0006 60.2281 17.5674 77.9675 31.827 70.4793 59.2555 50.7449 69.298 48.7498 62.7718 64.2383 52.4623 68.9526 52.915 54.6127 58.0381 29.4006 54.8979 26.1154 42.757 45.4043 8.14466 
57.5473 18.2616 111.201 59.6484 100.89 92.8962 49.7044 100.594 25.9789 77.3264 86.8 35.1278 122.447 44.8345 125.337 89.0185 108.606 113.885 100.021 114.737 108.577 102.463 111.181 94.1432 93.7586 93.7287 68.3643 90.327 70.2661 81.3108 93.586 10.6425 
54.1622 16.1773 107.098 53.7878 104.522 86.1288 66.08 98.0696 14.2164 83.3209 49.0063 46.8209 83.8659 22.8254 91.5101 59.206 78.128 87.6477 69.1584 95.139 85.1608 87.2046 104.407 80.5555 105.762 87.0567 85.7116 96.7177 55.6525 96.5171 44.7722 36.9821 
235.549 67.0697 473.97 228.994 485.525 382.402 345.356 467.619 115.008 451.122 148.12 333.498 346.407 157.888 440.062 137.287 407.776 292.704 270.074 385.848 112.994 379.541 199.868 286.421 334.678 174.58 379.202 197.955 317.932 301.266 180.423 146.708 
43.1605 15.8424 85.8862 43.6209 85.2092 69.1415 55.9523 80.8648 17.4742 73.6246 44.8917 51.239 78.8709 35.692 90.9621 55.1313 78.4146 79.6555 49.2306 90.4963 34.3458 86.2501 59.9026 76.2221 83.6899 74.7209 90.2009 82.9274 80.5402 87.4366 66.482 20.0367 
28.2095 7.71398 55.4914 25.9422 53.452 42.1514 34.1196 49.7661 23.2935 47.4382 44.6387 42.0641 62.0885 46.8741 62.6451 61.5197 47.1048 74.3249 32.8792 79.3493 49.9302 77.8338 75.8961 74.9692 90.6655 73.7162 89.8994 71.3042 78.4214 64.6463 68.2741 2.15828 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
];


figure,
bar3(input);view([90,-30,60]);set(gcf,'renderer','zbuffer'); set(gca, 'YLim', [0.0 32.0]);
set (gcf,'Position',[400,200,350,262], 'color','w')
% title('Spherical Harmonics Descriptors of The Origin Model'); 
xlabel('Frequency'); ylabel('Radius');zlabel('Energy');
x1=xlabel('Frequency');
set(x1,'Rotation',-60);

figure,
bar3(database_data);view([90,-30,60]);set(gcf,'renderer','zbuffer');  set(gca, 'YLim', [0.0 32.0]);
set (gcf,'Position',[400,200,350,262], 'color','w')
% title('Spherical Harmonics Descriptors of The Origin Model'); 
xlabel('Frequency'); ylabel('Radius');zlabel('Energy');
x1=xlabel('Frequency');
set(x1,'Rotation',-60);