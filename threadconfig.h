//
// Created by ylf9811 on 2021/3/20.
//

#ifndef YLF_BAM2FQ_THREADCONFIG_H
#define YLF_BAM2FQ_THREADCONFIG_H

class ThreadConfig {
public:
    ThreadConfig(int threadID);

public:
    int ThreadID;
    long TotalReads;
    long ReadsAligned;
    long ReadsQC;
    long TotalBases;
    long Bases[5];


};

#endif //YLF_BAM2FQ_THREADCONFIG_H
