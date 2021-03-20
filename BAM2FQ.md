# Bam QC

#### 0319



## base bam qc measures

- #### Total number of reads in the BAM file

  The `reads.total` QC metric counts the number of reads scanned in the BAM file(s). This includes reads later excluded due to low alignment score or as duplicates.
  The 42 BAMs contain 2.46 billion reads.

- #### Number of reads aligned to the reference genome

  The `reads.aligned` QC metric indicates the number of reads that were successfully aligned to the reference genome.
  The number of aligned reads is only 2% smaller, 2.42 billion.

- ##### Number of reads that passed minimum score filter

  At step 1 of RaMWAS, reads are filtered by the `scoretag` parameter, which is usually either the “MAPQ” field or “AS” tag in the BAM file. Reads with scores below `minscore` are excluded.

  The `reads.recorded` QC metric counts the number of reads that passed the score threshold.
  Almost of 2.2 billion reads passed the score threshold.


#### 0320

MAC

```
input : SRR_sort.bam

-o : 统计信息+输出
-b : 数行数
-  : 统计信息


-o
thread 1
➜  ylf_bam2fq time ./main $data/hg19/SRR_sort.bam tt.fq 1
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : tt.fq
thread : 1
benchmakr 0
is_write 1
total read is 6814517
total process is 6684853
total base is 668098535
total reads aligned is 6743512
A bases is 168454507
G bases is 163718457
C bases is 166439406
T bases is 169473729
N bases is 12436
./main $data/hg19/SRR_sort.bam tt.fq 1  4.81s user 0.86s system 112% cpu 5.024 total

thread 2
➜  ylf_bam2fq rm -rf tt.fq && time ./main $data/hg19/SRR_sort.bam tt.fq 2
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : tt.fq
thread : 2
benchmakr 0
is_write 1
total read is 6814517
total process is 6684853
total base is 668098535
total reads aligned is 6743512
A bases is 168454507
G bases is 163718457
C bases is 166439406
T bases is 169473729
N bases is 12436
./main $data/hg19/SRR_sort.bam tt.fq 2  4.84s user 0.77s system 217% cpu 2.576 total

thread 4
➜  ylf_bam2fq rm -rf tt.fq && time ./main $data/hg19/SRR_sort.bam tt.fq 4
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : tt.fq
thread : 4
benchmakr 0
is_write 1
total read is 6814517
total process is 6684853
total base is 668098535
total reads aligned is 6743512
A bases is 168454507
G bases is 163718457
C bases is 166439406
T bases is 169473729
N bases is 12436
./main $data/hg19/SRR_sort.bam tt.fq 4  5.57s user 0.85s system 408% cpu 1.570 total

-b
thread 1
➜  ylf_bam2fq rm -rf tt.fq && time ./main $data/hg19/SRR_sort.bam tt.fq 1 -b
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : tt.fq
thread : 1
benchmakr 1
is_write 0
total read is 6814517
totol process is 6684853
./main $data/hg19/SRR_sort.bam tt.fq 1 -b  2.47s user 0.19s system 106% cpu 2.493 total

thread 2
➜  ylf_bam2fq rm -rf tt.fq && time ./main $data/hg19/SRR_sort.bam tt.fq 2 -b
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : tt.fq
thread : 2
benchmakr 1
is_write 0
total read is 6814517
totol process is 6684853
./main $data/hg19/SRR_sort.bam tt.fq 2 -b  2.58s user 0.20s system 206% cpu 1.345 total

thread 4
➜  ylf_bam2fq rm -rf tt.fq && time ./main $data/hg19/SRR_sort.bam tt.fq 4 -b
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : tt.fq
thread : 4
benchmakr 1
is_write 0
total read is 6814517
totol process is 6684853
./main $data/hg19/SRR_sort.bam tt.fq 4 -b  2.53s user 0.18s system 403% cpu 0.671 total

-
thread 1
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : -
thread : 1
benchmakr 0
is_write 0
total read is 6814517
total process is 6684853
total base is 668098535
total reads aligned is 6743512
A bases is 168454507
G bases is 163718457
C bases is 166439406
T bases is 169473729
N bases is 12436
./main $data/hg19/SRR_sort.bam - 1  3.36s user 0.21s system 104% cpu 3.405 total

thread 2
➜  ylf_bam2fq rm -rf tt.fq && time ./main $data/hg19/SRR_sort.bam - 2
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : -
thread : 2
benchmakr 0
is_write 0
total read is 6814517
total process is 6684853
total base is 668098535
total reads aligned is 6743512
A bases is 168454507
G bases is 163718457
C bases is 166439406
T bases is 169473729
N bases is 12436
./main $data/hg19/SRR_sort.bam - 2  3.34s user 0.17s system 205% cpu 1.708 total

thread 4
➜  ylf_bam2fq rm -rf tt.fq && time ./main $data/hg19/SRR_sort.bam - 4
Starting Running
in : /Users/ylf9811/Desktop/QCQCQC/data/hg19/SRR_sort.bam
out : -
thread : 4
benchmakr 0
is_write 0
total read is 6814517
total process is 6684853
total base is 668098535
total reads aligned is 6743512
A bases is 168454507
G bases is 163718457
C bases is 166439406
T bases is 169473729
N bases is 12436
./main $data/hg19/SRR_sort.bam - 4  3.68s user 0.19s system 401% cpu 0.962 total
```

如果基本的信息统计，不输出文件的话，基本单线程慢1s。

6148

```
-
thread 1
ylf@gold6148:~/ylf_bam2fq$ rm -rf - && time ./main  /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam - 1
Starting Running
in : /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam
out : -
thread : 1
benchmakr 0
is_write 0
total read is 1212177023
total process is 1202406138
total base is 181563326838
total reads aligned is 1209374942
A bases is 53735158388
G bases is 36876338567
C bases is 37215422083
T bases is 53720088862
N bases is 16318938
real	11m38.142s
user	11m54.527s
sys	1m5.477s

thread 4
ylf@gold6148:~/ylf_bam2fq$ rm -rf - && time ./main  /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam - 4
Starting Running
in : /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam
out : -
thread : 4
benchmakr 0
is_write 0
total read is 1212177023
total process is 1202406138
total base is 181563326838
total reads aligned is 1209374942
A bases is 53735158388
G bases is 36876338567
C bases is 37215422083
T bases is 53720088862
N bases is 16318938
real	3m26.226s
user	13m52.322s
sys	0m30.249s

thread 8
ylf@gold6148:~/ylf_bam2fq$ rm -rf - && time ./main  /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam - 8
Starting Running
in : /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam
out : -
thread : 8
benchmakr 0
is_write 0
total read is 1212177023
total process is 1202406138
total base is 181563326838
total reads aligned is 1209374942
A bases is 53735158388
G bases is 36876338567
C bases is 37215422083
T bases is 53720088862
N bases is 16318938
real	1m44.212s
user	14m0.762s
sys	0m27.920s

thread 16
ylf@gold6148:~/ylf_bam2fq$ rm -rf - && time ./main  /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam - 16
Starting Running
in : /home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam
out : -
thread : 16
benchmakr 0
is_write 0
total read is 1212177023
total process is 1202406138
total base is 181563326838
total reads aligned is 1209374942
A bases is 53735158388
G bases is 36876338567
C bases is 37215422083
T bases is 53720088862
N bases is 16318938
real	1m0.362s
user	16m11.765s
sys	0m28.236s
```

