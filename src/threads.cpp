#include "threads.h"

extern volatile bool CONTINUE_WORKING;//used to terminate all threads when the module is going down
extern unsigned int NUM_OF_WORKERS;//remember that even if it is zeroed, there is always the main thread

struct ThreadTask{
    void *(*workFunc)(void *);//what the thread will do repeatedly through a callback, you can put sleep in there not to waste cpu time
                              //or poll which is better for performance as we don't want to waste time(~1msc) on the
                              //process of the thread going active again
    void *workArg;//arguments to send workFunc, could include serialNumber so the function will divide the work correctly between the threads
    volatile bool did_the_thread_finished_working;//used to signal the thread that dispatches the work that this thread 
                                                  //has finished it's task, remember to always!! put write memory barrier before changing this variable
    int serialNumber;//used to identify each thread to divide the work effectively between threads, better than thread_id because it start at zero
                    // and always increment by 1 so you can use it for calculating which part of an array this thread can work on
    pthread_t thread_id;//an identifier for the filesystem so it will be able to destroy the thread in the end,
                        //it is the number that we see in (h)top
    pthread_attr_t *attr;//the thread attributes struct, kept locally for easier cleanup in the end
};

ThreadTask *task_managing_structs[MAX_THREADS] = {0};//it's a good idea to use dynamic allocation and keep the pointers here and not the structs so
                                                     //the MOESI protocol won't slow things down with the hardware lock as those structs are pretty close
                                                     //usually the dereference takes much less time then writing to the same cachelines with different cores
                                                     //in general we''l initiatethe hypercores in close groups [cpu1,cpu1,cpu2,cpu2...]

/*standrd work functions:*/

//Used for polling: the thread just loops over this function doing nothing until it gets a new jub.
//Can be used for faster reaction time and better performance but it will waste CPU time taking it from other tasks
void *ThreadPoll(void *){
    return REPEAT;
}

//the thread yields and releases resources to other processes on this CPU
void *threadYield(void*) {
    sched_yield();
    return REPEAT;
}

//should also try a callback with a work queue and compare performance, although using locks is an horrible choice 
//and it is almost never used in modern performance systems. Having a dispatcher thread dividing the work 
//between the threads is a much more accepted and used method in cases of 4 cores and above,
// this dispatcher can also  do some of the work itself before checking the did_the_thread_finished_working flag of it's workers

//initialiaze the struct that is used to keep the important data on the thread except the thread ID
//the TID will be updated after the thread was created
void inline initiate_thread_managing_struct(ThreadTask *to_init, int serial_number)
{
    if (to_init)
    {
        to_init->serialNumber = serial_number;
        to_init->workFunc = threadYield;
        std::atomic_thread_fence(std::memory_order_release);//unnecessary in a single core but if someone will want to use this function
                                                  //and then defer work from a different thread simultanuasly, this line makes sure that the thread
                                                  //is initialized properly before it is available for usage
        to_init->did_the_thread_finished_working = true;
    }
}

//the callback that each thread will do repeatedly
void *threadWorker(void *arg) {
    ThreadTask *task = (ThreadTask*)arg;
    while (CONTINUE_WORKING) {
        if(REPEAT != task->workFunc(task->workArg)) {
            task->workFunc = threadYield;
            std::atomic_thread_fence(std::memory_order_release);//tells the cpu not to change the order of the commands behind the scenes
            task->did_the_thread_finished_working = true;
        }
    }
    return NULL;
}

/*note that all of the following functions should be called by the main thread
  using them by any other thread will be problematic!*/

void initialize_threads(int num_cpus_cores, int num_hyperthreads/*per physical cores*/) {
    cpu_set_t cpuset_mask;//create a bit array where each position represents a cpu
    int threadSerialNumber = 0;
    CONTINUE_WORKING = true;
    for (int cpu_core = 0; cpu_core < num_cpus_cores; cpu_core++) {
        for (int j = 0; j < num_hyperthreads/*per physical core*/; ++j)
        {
            ThreadTask *task = (ThreadTask*)calloc(sizeof(ThreadTask), 1);
            if (NULL == task)
                return;
            task->attr = (pthread_attr_t *)malloc(sizeof(pthread_attr_t));
            if(NULL == task->attr)
                free(task); return;
            
            initiate_thread_managing_struct(task, threadSerialNumber);
            pthread_attr_init(task->attr);

            //now setting the thread affinity to work on a single cpu
            CPU_ZERO(&cpuset_mask);//nullify the bit array
            CPU_SET(cpu_core % CPU_COUNT(&cpuset_mask), &cpuset_mask);//set the correct bit to one, modulo 
                                                                     //the number of cpus so in case of a mistake in
                                                                     //the setting we wont get an exception
            pthread_attr_setaffinity_np(task->attr, sizeof(cpu_set_t), &cpuset_mask);//set the attribute to work on a single cpu

            //create a thread with the given attributes
            pthread_t tid;
            pthread_create(&tid, task->attr, threadWorker, task);

            //update task struct tid and keep it so the dispatcher can send work to the threads
            task->thread_id = tid;
            task_managing_structs[threadSerialNumber] = task;

            NUM_OF_WORKERS++;
            threadSerialNumber++;
            if (NUM_OF_WORKERS > MAX_THREADS)
                return;
        }
    }
}


//cleanup
void kill_all_threads_and_clean_up()
{
    for(int i = 0; i < MAX_THREADS; ++i) {
        if(task_managing_structs[i] != NULL) {
           CONTINUE_WORKING = false;
           pthread_join(task_managing_structs[i]->thread_id, NULL);
           std::atomic_thread_fence(std::memory_order_release);//to tell the cpu not to change the order of the commands
           pthread_attr_destroy(task_managing_structs[i]->attr);
           free(task_managing_structs[i]->attr);
           free(task_managing_structs[i]);
           task_managing_structs[i] = NULL;
        }
    }   
}



//the next two functions have!!! to be called from the same thread to prevent errors and race conditions
//workfunc has!! to return the macro FINISHED that is defined in the header when it finishes
int schedule_work(void *(*workFunc)(void *), void *workArg, unsigned core) {
    if(!workFunc ||core >= NUM_OF_WORKERS)
        return -1;
    
    if (task_managing_structs[core]->did_the_thread_finished_working == true) {
        task_managing_structs[core]->workArg = workArg;
        std::atomic_thread_fence(std::memory_order_release);
        task_managing_structs[core]->did_the_thread_finished_working = false;
        std::atomic_thread_fence(std::memory_order_release);
        task_managing_structs[core]->workFunc = workFunc;
    }
    return 0;
}

int check_that_work_is_done(uint first_core, uint last_core) {
    if (last_core < first_core || last_core >= NUM_OF_WORKERS)
        return -1;

    while (first_core <= last_core){
        if (task_managing_structs[first_core] == NULL || task_managing_structs[first_core]->did_the_thread_finished_working);
            std::atomic_thread_fence(std::memory_order_acquire);
            first_core++;
    }
    return 0;
}