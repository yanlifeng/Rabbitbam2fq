

## BWA的使用

DOTO

- [x] mem和sw算法的输出结果不一样
- [x] aln输出二进制文件？
- [x] map文件是啥
- [x] sam中第二列是0？
- [x] main.c改为cpp
- [ ] 加上fprintf
- [x] htslib icc -O3
- [x] seq和qual转化时压位 单线程还能快0.5？
- [ ] CRC32 checksum ？？
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

0221 optimize nibble2basere function, now single is 0.5s faster.

| threadNumber | Cost                                         |
| ------------ | -------------------------------------------- |
| 1            | 4.97s user 0.94s system 132% cpu 4.449 total |
| 2            | 4.95s user 0.84s system 252% cpu 2.287 total |
| 4            | 5.82s user 1.00s system 498% cpu 1.369 total |
| 8            | 6.51s user 1.08s system 632% cpu 1.199 total |

0222 摸大🐟

0223 it seems that read data from sam file is faster.(maybe the reason is that read part can be muti-threading)

| threadNumber | Cost                                         |
| ------------ | -------------------------------------------- |
| 1            | 2.65s user 1.11s system 160% cpu 2.343 total |
| 2            | 2.67s user 1.11s system 296% cpu 1.276 total |
| 4            | 4.10s user 1.54s system 550% cpu 1.025 total |
| 8            | 3.96s user 1.42s system 540% cpu 0.996 total |

今天又有新发现，bam_read1里面也有多线程。在bgzf_read_block中使用了多线程，具体的还在看。



现在的思路有点迷糊了，view的多线程看上去效果很好，具体怎么实现的（for unidex bam file）需要好好看看。会不会它实际上也是外层串行的？对于index的数据还能更快？

（https://www.cnblogs.com/hanyonglu/archive/2011/05/07/2039916.html学了新的传参方法）

通过分析文件名，得到format.format是bam；通过分析前几个字节的信息，得到format.compression是bgzf，然后调用bgzf_thread_pool函数弄了一个线程池。

```c
pthread_create(&mt->io_task, NULL,fp->is_write ? bgzf_mt_writer : bgzf_mt_reader, fp);
```

这个是创建线程的语句

每次从任务队列（hts_tpool_process q）里面拿出一个result

```
/*
 * An output, after job has executed.
 */
struct hts_tpool_result {
    struct hts_tpool_result *next;
    void (*result_cleanup)(void *data);
    uint64_t serial; // sequential number for ordering
    void *data;      // result itself
};
```

res中的data，即bgzf_job，然后就有了想要的指针什么的

```
typedef struct bgzf_job {
    BGZF *fp;
    unsigned char comp_data[BGZF_MAX_BLOCK_SIZE];
    size_t comp_len;
    unsigned char uncomp_data[BGZF_MAX_BLOCK_SIZE];
    size_t uncomp_len;
    int errcode;
    int64_t block_address;
    int hit_eof;
} bgzf_job;
```

现在的关键是找到hts_tpool_process q这个队列是怎么构建的。

不太对，最最开始的hts_tpool_init里面

```
typedef struct {
    struct hts_tpool *p;
    int idx;
    pthread_t tid;
    pthread_cond_t  pending_c; // when waiting for a job
} hts_tpool_worker;
```

```
for (t_idx = 0; t_idx < n; t_idx++) {
    printf("new thread %d\n", t_idx);
    hts_tpool_worker *w = &p->t[t_idx];
    p->t_stack[t_idx] = 0;
    w->p = p;
    w->idx = t_idx;
    pthread_cond_init(&w->pending_c, NULL);
    if (0 != pthread_create(&w->tid, NULL, tpool_worker, w)) {
        goto cleanup;
    }
}
```

这个地方，在啥都没干的时候定义了n个线程，任务就是static void *tpool_worker(void *arg)；这个时候所有的线程就已经开始干活了，

（https://www.cnblogs.com/qyaizs/articles/2039101.html关于tpyedef struct等）

实际上这个时候已经有n个线程开始执行tpool_worker了，我们去tpool_worker里面看一下。

```
* Once woken, each thread checks each process-queue in the pool in turn,
* looking for input jobs that also have room for the output (if it requires
* storing).  If found, we execute it and repeat.
```

这个是注释，基本上也说清楚了，就是检测process-queue中是否有可以操作的job，有就进行hts_tpool_add_result，不过同时运行了j->func(j->arg)，暂时不知道这个func是啥，不过应该是热点，因为它的返回值就是data，后面对data的处理就是上面samloop中bgzf_mt_reader的简单拷贝了。

然后就和上面发现的衔接起来了，hts_set_opt中调用了bgzf_thread_pool，创建了io_task线程，运行bgzf_mt_reader方法：

```
pthread_create(&mt->io_task, NULL, fp->is_write ? bgzf_mt_writer : bgzf_mt_reader, fp);
```

里面的热点，即vtune里面飙出来的bgzf_decode_func：

```
if (hts_tpool_dispatch3(mt->pool, mt->out_queue, bgzf_decode_func, j,
                        job_cleanup, job_cleanup, 0) < 0) {
    job_cleanup(j);
    goto err;
}
```

hts_tpool_dispatch3：Adds an item to the work pool.

```
j->func = exec_func;
j->arg = arg;
j->job_cleanup = job_cleanup;
j->result_cleanup = result_cleanup;
j->next = NULL;
j->p = p;
j->q = q;
j->serial = q->curr_serial++;
```

```
if (q->input_tail) {
    q->input_tail->next = j;
    q->input_tail = j;
} else {
    q->input_head = q->input_tail = j;
}
```

可以看到func即bgzf_decode_func，现在把j放到了q里面，同时tpool_worker检测到j，执行j->func，我们回到tpool_worker里面：

```
pthread_mutex_unlock(&p->pool_m);

DBG_OUT(stderr, "%d: Processing queue %p, serial %"PRId64"\n",
        worker_id(j->p), q, j->serial);

if (hts_tpool_add_result(j, j->func(j->arg)) < 0)
    goto err;
//memset(j, 0xbb, sizeof(*j));
free(j);

pthread_mutex_lock(&p->pool_m);
```

可以看到j->func(j->arg)是没有锁的，即多线程执行的！



0224

今天想着去ASC上跑一跑，结果报错了

```
./fast: symbol lookup error: ./fast: undefined symbol: fq_write1
```

搜了一下，大体的意思就是链接库的版本不对，我新编译好的带fq_write1函数的版本没用上，ldd(otool -L on mac)命令发现

```
	libhts.so.3 => /usr/local/lib/libhts.so.3 (0x00007fcbe9a1f000)
```

不知道为啥不去我指定的目录里面找，还是去了默认的目录，暂时的解决办法是往默认目录重新编译一遍。

asc上跑的很慢很慢的，最后好像要释放一些东西还是咋，大概是hdd读写太慢了。

0225

分析vtune发现

单线程的

![image-20210225102501724](/Users/ylf9811/Library/Application Support/typora-user-images/image-20210225102501724.png)

不开线程的

![image-20210225102540234](/Users/ylf9811/Library/Application Support/typora-user-images/image-20210225102540234.png)

可以看到，不开-@的比较正常，读写三七开，读里面bgzf_read函数里面的read_block()最慢，和预期的一样；

但是开-@的程序，哪怕是只开一个线程，时间都记在了thread里面。