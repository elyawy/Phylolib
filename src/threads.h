#ifndef __threads
#define __threads

#include <pthread.h>
#include <atomic>
#include <stdlib.h>
#include <unistd.h>

/*
 * This file implements a dispatcher consumer thread module. The idea is to let the program have a main thread that
 * do most of the work and when a profiler is detecting that some function is slowing the code down,
 * you can use this module to throw more cores on it without a lot of complications.
 */

#define MAX_THREADS 40 //2*max num of cores, I suggest keeping the number of threads 2 off this max to keep a physical core for the executing core 
#define NUM_HYPERTHREADS 1//useful in case you want the cores to work on the same areas and has the same cache
unsigned int NUM_OF_WORKERS = 0;//remember that even if it is zeroed, there is always the main thread
#define REPEAT (void *)0xDEADBEEF //used to signal this model that the work function should be called again, any other return value
                                  //will be interpreted as exit

volatile bool CONTINUE_WORKING = true;//used to terminate all threads when the module is down
typedef struct ThreadTask ThreadTask;

extern ThreadTask *task_managing_structs[MAX_THREADS];//it's a good idea to use dynamic allocation and keep the pointers here and not the structs so
                                                //the MOESI protocol won't slow things down with the hardware lock as those structs are pretty close
                                              //usually the dereference takes much less time then writing to the same cachelines with different cores
                                              //in general we''l initiatethe hypercores in close groups [cpu1,cpu1,cpu2,cpu2...]


/*note that all of the following functions should be called by the main thread
  using them by any other thread will be problematic!*/

void initialize_threads(int num_cpus_cores, int num_hyperthreads/*per physical cores*/);
//cleanup
void kill_all_threads_and_clean_up();

//the next two functions have!!! to be called from the same thread to prevent errors and race conditions
int schedule_work(void *(*workFunc)(void *), void *workArg, uint first_core, uint last_core);
int check_that_work_is_done(uint first_core, uint last_core);

#endif 