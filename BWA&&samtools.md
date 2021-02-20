

## BWA的使用

DOTO

- [x] mem和sw算法的输出结果不一样
- [x] aln输出二进制文件？
- [x] map文件是啥
- [x] sam中第二列是0？
- [x] main.c改为cpp
- [ ] 加上fprintf
- [ ] htslib icc -O3

- [ ] seq和qual转化时压位
- [ ] 
- [ ] 

http://bio-bwa.sourceforge.net/bwa.shtml

### 安装

```c++
git clone https://github.com/lh3/bwa.git
cd bwa
make
```

### 建立索引

```c++
./bwa index -a bwtsw hg19.fa（./bwa index hg19.fa）～3000s
生成这几个文件(6148. ~/data/hg19/)
hg19.fa  hg19.fa.amb  hg19.fa.ann  hg19.fa.bwt  hg19.fa.pac  hg19.fa.sa
```

### 比对

#### BWA-backtrack

```
./bwa aln hg19.fa test.fq > alnTest.sam
```

（test.fa长度大于100，可能会出现一丢丢问题）

～40s

#### BWA-SW

```c++
./bwa bwasw hg19.fa test.fq > test.sam
```

～20s

#### BWA-MEM

```c++
./bwa mem hg19.fa test.fq > memTest.sam
```

～10s

### sam文件解析

生成的对比信息都存在sam文件中(bam是sam的二进制形式)。

##### 头部注释部分

- @HD VN:1.0 SO:unsorted （排序类型）
  头部区第一行：VN是格式版本；SO表示比对排序的类型，有unknown（default），unsorted，queryname和coordinate几种。samtools软件在进行行排序后不能自动更新bam文件的SO值，而picard却可以。
- @SQ SN:contig1 LN:9401 （序列ID及长度）
  参考序列名，这些参考序列决定了比对结果sort的顺序，SN是参考序列名；LN是参考序列长度；每个参考序列为一行。
  例如：@SQ SN:NC_000067.6 LN:195471971
- @RG ID:sample01 （样品基本信息）
  Read Group。1个sample的测序结果为1个Read Group；该sample可以有多个library的测序结果，可以利用bwa mem -R 加上去这些信息。
  例如：@RG ID:ZX1_ID SM:ZX1 LB:PE400 PU:Illumina PL:Miseq
  ID：样品的ID号 SM：样品名 LB：文库名 PU：测序以 PL：测序平台
  这些信息可以在形成sam文件时加入，ID是必须要有的后面是否添加看分析要求
- @PG ID:bowtie2 PN:bowtie2 VN:2.0.0-beta7 （比对所使用的软件及版本）
  例如：@PG ID:bwa PN:bwa VN:0.7.12-r1039 CL:bwa sampe -a 400 -f ZX1.sam -r @RG ID:ZX1_ID SM:ZX1 LB:PE400 PU:Illumina PL:Miseq …/0_Reference/Reference_Sequence.fa ZX_HQ_clean_R1.fq.sai ZX_HQ_clean_R2.fq.sai …/2_HQData/ZX_HQ_clean_R1.fq …/2_HQData/ZX_HQ_clean_R2.fq
  这里的ID是bwa，PN是bwa，VN是0.7.12-r1039版本。CL可以认为是运行程序@RG是上面RG表示的内容，后面是程序内容，这里的@GR内容是可以自己在运行程序是加入的

##### 比对结果部分

比对结果部分每行标示一个read与参考序列的比对信息，前11列为必须字段，顺序固定，其余列是可选字段。

- 第一列Query Name：read的名称，即片段的编号。

- 第二列FLAG：如果不是以下数字中的一个，则是一下数据某几个的和：
  ```bash
  - 1：标示对应的二进制为01，标示read有多个测序数据，一般理解为有双端测序数据，另一条没有过滤掉；
  - 2：二进制为10，标示read的多个片段都有比对结果，双端的read都比对上了；
  - 4：表示这条read没有比对上；
  - 8：标示下一条read没有比对上；
  - 16：表示这条read的反向比对上了；
  - 32：表示这条read的下一条的反向没有比对上；
  - 64：表示样本中第一个片段；
  - 128：表示样本中最后一条片段；
  - 256：表示第二次比对；
  - 512：表示比对的质量不合格；
  - 1204：表示read是pcr或光学副本产生的；
  - 2048：表示辅助比对结果；
  ```

- 第三列Reference Name：参考序列的名称，或者比对到参考序列上的染色体号。比对不上为*。

- 第四列Position：比对上的位置，从1开始计数（顺着链的方向从1数起，哪个位置开始匹配），没有比对上为0。

- 第五列Mapping Quality：比对的质量分数，越高表示比对的越准确。

- 第六列CIGAR：表示比对的结果。

  ```rust
  * M：表示match或mismatch
  * I：表示插入
  * D：表示删除
  * N：表示skipped，跳过这段区域
  * S：表示被剪切的序列存在于序列中
  * H：表示被剪切的序列不存在于序列中
  * P：表示padding（填补）
  * =：表示match
  * X：表示mismatch（错配，位置是一一对应的）
  ```

- 第七列RNEXT：表示下一个片段比对上的参考序列的编号，比对不上用’*‘，该片段和下一个片段比对上同一个参考片段，用’=‘。

- 第八列PNEXT：表示下一个片段比对上的位置，如果不可用，此处为0。

- 第九列TLEN：表示Template的长度。如果第八列大于第四列，则为正数，否则负数。

- 第十列SEQ：表示序列片段的序列信息，（注意CIGAR中M/I/S/=/X对应数字的和要等于序列长度），表示read的碱基序列，如果是比对到互补链上则是反转互补序列。

- 第十一列QUAL：表示read的质量，用ASCII编码表示。

### 验证

经过自己手动验证，sam文件中的染色体编号和序列位置找的基本是正确的。



## Samtools的使用

### 安装

```
下载 http://www.htslib.org/download/
cd samtools-1.x    # and similarly for bcftools and htslib
./configure --prefix=/where/to/install
make
make install
export PATH=/where/to/install/bin:$PATH    # for sh or bash users
```

### 使用

#### view

```
samtools view [options] in.bam|in.sam|in.cram [region...]
```

#### index

```
samtools index [-bc] [-m INT] aln.bam|aln.cram
```

#### sort

```
samtools sort [-l level] [-m maxMem] [-o out.bam] [-O format] [-n] -T out.prefix [-@ threads] [in.bam]
```

bam2fq

```
samtools bam2fq [-nO] [-s <outSE.fq>] <in.bam>
```



(新CPU的TDP按照150算，可能算小了，没事到时候多出来的节点可以不用)

NF5280M6：150x2+7.5x16+10+10=440w

| Mode     | Name                  | Power consumption                        |
| -------- | --------------------- | ---------------------------------------- |
| CPU Mode | NF5280M6              | 440w*7                                   |
|          | NVIDIA V100s-PCIE-32G | 20w*4（待机功率，不确定）                |
|          | GbE switch            | 30w                                      |
|          | HDR-IB switch         | 130w（这个switch不知道是都用还是选一个） |
|          | Total                 | 3320                                     |

| Mode     | Name                  | Power consumption                        |
| -------- | --------------------- | ---------------------------------------- |
| GPU Mode | NF5280M6              | 160w*7（待机功率，不确定）               |
|          | NVIDIA V100s-PCIE-32G | 250*4                                    |
|          | GbE switch            | 30w                                      |
|          | HDR-IB switch         | 130w（这个switch不知道是都用还是选一个） |
|          | Total                 | 2280w                                    |



### RabbitBam2Fq

test.fq：557M 2000000reads

|                              | time(Mac) | Time(6148) |
| ---------------------------- | --------- | ---------- |
| samtools                     | 3.380     |            |
| htslib printf >              | 15.756    |            |
| htslib no output             | 2.058     |            |
| htslib ofstream(no buffer)   | 4.273     |            |
| htslib ofstream(100M buffer) | 2.225     |            |
| * -O3 -ffast-math            | 1.730     |            |
|                              |           |            |

Index:

|                              | time(Mac) | Time(6148) |
| ---------------------------- | --------- | ---------- |
| samtools                     | 3.387     |            |
| htslib printf >              | x         |            |
| htslib no output             | x         |            |
| htslib ofstream(no buffer)   | x         |            |
| htslib ofstream(100M buffer) | 2.437     |            |

SRR_sort.fq : 1.4G SRR_sort.bam 300M  6684853 reads -O3 -ffast-math

|                                | Time(MAC)       |
| ------------------------------ | --------------- |
| Samtools                       | 11.415          |
| Htslib just read and parser    | 3.792           |
| Htslib ofstream(100M buffer)   | 5.431           |
| Htslib just read but no parser | 3.553 ????????? |
|                                |                 |



0218版本main.c SRR_sort.bam 313MB

| threadNumber | Cost                                          |
| ------------ | --------------------------------------------- |
| 1            | 5.54s user 2.74s system 97% cpu 8.535 total   |
| 2            | 5.56s user 2.46s system 165% cpu 4.835 total  |
| 4            | 5.77s user 2.39s system 280% cpu 2.914 total. |

0219版本test.cpp ... 优化了输出 保证了输出的正确性（仍然乱序）

| threadNumber | Cost                                         |
| ------------ | -------------------------------------------- |
| 1            | 4.89s user 0.78s system 98% cpu 5.769 total  |
| 2            | 5.15s user 0.81s system 167% cpu 3.566 total |
| 4            | 5.45s user 0.85s system 296% cpu 2.123 total |
| 8            | 7.11s user 1.13s system 507% cpu 1.623 total |

现在尝试优化单线程的版本：

| version       | cost                                        |
| ------------- | ------------------------------------------- |
| init main.cpp | 5.18s user 0.62s system 97% cpu 5.941 total |
|               |                                             |
|               |                                             |
|               |                                             |

samtool view ！！

| threadNumber | Cost                                         |
| ------------ | -------------------------------------------- |
| 1            | 5.85s user 1.23s system 135% cpu 5.241 total |
| 2            | 5.84s user 1.08s system 254% cpu 2.716 total |
| 4            | 6.97s user 1.16s system 485% cpu 1.674 total |
| 8            | 7.84s user 1.06s system 649% cpu 1.371 total |

samtool view ！！ no write

| threadNumber | Cost                                         |
| ------------ | -------------------------------------------- |
| 1            | 4.40s user 0.43s system 121% cpu 3.961 total |
| 2            | 4.56s user 0.46s system 244% cpu 2.049 total |
| 4            | 4.96s user 0.33s system 479% cpu 1.103 total |
| 8            | 5.66s user 0.27s system 675% cpu 0.877 total |

0220 new bam2fq (almost use htslib test_view.c)

| threadNumber | Cost                                         |
| ------------ | -------------------------------------------- |
| 1            | 4.99s user 0.91s system 133% cpu 4.422 total |
| 2            | 5.08s user 0.91s system 257% cpu 2.327 total |
| 4            | 5.86s user 1.00s system 505% cpu 1.356 total |
| 8            | 6.47s user 1.09s system 608% cpu 1.241 total |

0220 new bam2fq (almost use htslib test_view.c)  and  get right fq file（include reverse some reads）

| threadNumber | Cost                                         |
| ------------ | -------------------------------------------- |
| 1            | 5.61s user 1.09s system 132% cpu 5.052 total |
| 2            | 5.54s user 0.93s system 251% cpu 2.579 total |
| 4            | 6.52s user 1.05s system 485% cpu 1.560 total |
| 8            | 7.28s user 1.00s system 628% cpu 1.316 total |