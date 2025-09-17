I was unable to increase the stack size. as it results in the following error.

ulimit: value exceeds hard limi

Instead, in main.py, I increated sys.setrecursionlimit(1000000) sys.setrecursionlimit(10**7) to to handle data4.

This allowed my device to run through data4 with no problems.


(base) spencerlockwood@w134-87-061-197 genome-assembly % python main.py data1
short_1.fasta 8500 100
short_2.fasta 8500 100
long.fasta 250 1000
0 15650
1 9997
2 9997
3 9990
4 9990
5 9956
6 9956
7 4615
8 3277
9 828
10 684
11 669
12 669
13 666
14 666
15 655
16 654
17 639
18 639
19 636

(base) spencerlockwood@w134-87-061-197 genome-assembly % python main.py data2
short_1.fasta 5000 100
short_2.fasta 5000 100
long.fasta 500 1000
0 15744
1 10013
2 10013
3 9992
4 9992
5 9992
6 5752
7 5171
8 4664
9 3309
10 1009
11 938
12 829
13 733
14 654
15 652
16 652
17 652
18 652
19 652

(base) spencerlockwood@w134-87-061-197 genome-assembly % python main.py data3
short_1.fasta 2500 100
short_2.fasta 2500 100
long.fasta 500 1000
0 9824
1 9824
2 9824
3 9824
4 9824
5 9824
6 9824
7 9824
8 3656
9 3656
10 3592
11 3592
12 2604
13 1848
14 1654
15 1517
16 1431
17 1408
18 1352
19 1239

(csc427_venv) (base) spencerlockwood@w134-87-061-197 genome-assembly % python main.py data4                     
short_1.fasta 25000 100
short_2.fasta 25000 100
long.fasta 5000 1000
0 173867
1 173801
2 159255
3 152869
4 25669
5 21727
6 21727
7 18827
8 18827
9 10981
10 6798
11 6798
12 6467
13 6423
14 5715
15 5715
16 5219
17 5219
18 4798
19 4798


I am getting the following:

For data1: N50 = 9990
For data2: N50 = 9992
for data3:  N50 = 9824