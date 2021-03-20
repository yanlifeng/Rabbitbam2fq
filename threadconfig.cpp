//
// Created by ylf9811 on 2021/3/20.
//
#include "threadconfig.h"

ThreadConfig::ThreadConfig(int threadID) {
    ThreadID = threadID;
    TotalReads = 0;
    ReadsAligned = 0;
    ReadsQC = 0;
    TotalBases = 0;
    for (int i = 0; i < 5; i++)
        Bases[i] = 0;
};
