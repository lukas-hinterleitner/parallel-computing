#include <iostream>
#include <thread>
#include <chrono>
#include <fstream>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <future>
#include "a1-helpers.hpp"

template <typename T>
class SafeQ
{
private:
    std::queue<T> q; // no other data structures are allowed
    std::mutex m;
    std::condition_variable cv;
public:
    void push(T value)
    {
        m.lock();
        q.push(value);
        m.unlock();
        cv.notify_one();
    }

    void pop(T &value)
    {
        // if not empty, remove the element from the queue
        m.lock();
        if (!q.empty())
        {
            value = q.front();
            q.pop();
        }
        m.unlock();
    }

    T wait_and_pop()
    {
        std::unique_lock lock(m);
        cv.wait(lock, [this]{ return !q.empty(); });
        const int value = q.front();
        q.pop();
        return value;
    }

    size_t size()
    {
        std::unique_lock lock(m);
        return q.size();
    }

    bool empty()
    {
        std::unique_lock lock(m);
        return q.empty();
    }
};

bool producer_finished = false;
std::mutex producer_finished_mutex;

/**
 * To be executed by the master thread
 * The function reads numbers from the file
 * and puts them into the given queue
 *
 * @param[in] filename
 * @param[inout] q
 * @returns The number of produced items
 *
*/
int producer(const std::string& filename, SafeQ<int> &q)
{
    int produced_count = 0;

    // while there are entries in the file
    // put numbers into the queue "q"
    std::ifstream ifs(filename);

    while (!ifs.eof()) {
        int num;
        ifs >> num;
        q.push(num);
        produced_count++;
    }

    producer_finished_mutex.lock();
    producer_finished = true;
    producer_finished_mutex.unlock();

    ifs.close();

    return produced_count;
}

/**
 * To be executed by worker threads
 * The function removes a number from the queue "q"
 * and does the processing
 * Implement 2 versions with atomic and with mutexes
 * extend as needed
 *
 * @param[inout] q
 * @param[inout] primes
 * @param[inout] nonprimes
 * @param[inout] mean
 * @param[inout] number_counts
 *
*/
void worker(SafeQ<int> &q, int &primes, int &nonprimes, double &sum, int &consumed_count, std::vector<int> &number_counts, std::mutex &metrics_mutex)
{
    // implement: use synchronization
    // Note: This part may need some rearranging and rewriting
    // the while loop cannot just check for the size,
    // it has to now wait until the next element can be popped,
    // or it has to terminate if producer has finished and the queue is empty.

    for(;;)
    {
        int num;
        {
            producer_finished_mutex.lock();
            if (producer_finished && q.empty())
            {
                producer_finished_mutex.unlock();
                break;
            }

            num = q.wait_and_pop();
            producer_finished_mutex.unlock();
        }

        const bool is_prime = kernel(num) == 1;

        metrics_mutex.lock();
        consumed_count++;
        number_counts[num % 10]++;
        sum += num;

        if (is_prime) primes++;
        else nonprimes++;
        metrics_mutex.unlock();
    }
}

int main(int argc, char **argv)
{
    int num_threads = (int) std::thread::hardware_concurrency(); // you can change this default to thread::hardware_concurrency()
    bool no_exec_times = false, only_exec_times = false;; // reporting of time measurements
    std::string filename = "input.txt";
    parse_args(argc, argv, num_threads, filename, no_exec_times, only_exec_times);

    // The actual code
    int primes, nonprimes = 0;
    int consumed_count = 0;
    double mean, sum = 0.0;

    // vector for storing numbers ending with different digits (0-9)
    std::vector<int> number_counts(10, 0);

    std::mutex metrics_mutex;

    // Queue that needs to be made safe
    // In the simple form it takes integers
    SafeQ<int> q;

    // put you worker threads here

    std::vector<std::thread> workers;

    // time measurement
    auto t1 = std::chrono::high_resolution_clock::now();

    // implement: call the producer function with futures/async
    std::future<int> produced_count_future = async(std::launch::async, producer, filename, std::ref(q));

    // implement: spawn worker threads - transform to spawn num_threads threads and store in the "workers" vector
    for (int i=0;i<num_threads;++i) {
        workers.emplace_back(worker, std::ref(q), std::ref(primes), std::ref(nonprimes), std::ref(sum),
                             std::ref(consumed_count), ref(number_counts), ref(metrics_mutex));
    }

    for (int i = 0; i < num_threads; ++i) {
        workers[i].join();
    }

    mean = sum/consumed_count;
    // end time measurement
    auto t2 = std::chrono::high_resolution_clock::now();

    int produced_count = produced_count_future.get();

    // do not remove
    if ( produced_count != consumed_count ) {
        std::cout << "[error]: produced_count (" << produced_count << ") != consumed_count (" << consumed_count << ")." << std::endl;
    }

    // priting the results
    print_output(num_threads, primes, nonprimes, mean, number_counts, t1, t2, only_exec_times, no_exec_times);

    return 0;
}