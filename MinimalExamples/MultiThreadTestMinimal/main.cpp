#include <iostream>
#include <thread>

//This function will be called from a thread

void call_from_thread()
{
    std::cout << "Hello, World" << std::endl;
}

int main()
{
    //Launch a thread
    int k = 10000;
    std::thread t[k];
    for (int i = 0; i < k; i++)
    {
        t[i] = std::thread(call_from_thread);
    }
    //Join the thread with the main thread
    for (int i = 0; i < k; i++)
    {
        t[i].join();
    }
   // t1.join();

    return 0;
}
