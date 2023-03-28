#include <iostream>
#include <thread>
#include <chrono>
#include <fstream>
#include <queue>
#include <condition_variable>
#include <future>
#include <atomic>
#include "a1-helpers.hpp"

template <typename T>
class SafeQ
{
private:
    std::queue<T> q; // no other data structures are allowed
    std::mutex m;
    std::condition_variable cv;
    // extend as needed
public:
    void push(T value)
    {
        // synchronization needed?
        m.lock();
        q.push(value);
        m.unlock();
        cv.notify_one();

    }
    void pop(T &value)
    {
        // todo:
        // in a thread-safe way take the front element
        // and pop it from the queue
        // multiple consumers may be accessing this method

        // if not empty, remove the element from the queue
        std::unique_lock<std::mutex> lock(m);
        if (!q.empty()) {
            value = q.front();
            q.pop();
        }
    }

    T wait_and_pop()
    {
        // todo:
        // in a thread-safe way take the front element
        // and pop it from the queue
        // multiple consumers may be accessing this method

        std::unique_lock<std::mutex> lock(m);
        cv.wait(lock, [this]{ return !q.empty(); });
        const int value = q.front();
        q.pop();
        return value;
    }

    size_t size()
    {
        std::unique_lock<std::mutex> lock(m, std::defer_lock);
        return q.size();
    }

    bool empty()
    {
        std::unique_lock<std::mutex> lock(m, std::defer_lock);
        return q.empty();
    }
};

std::atomic<bool> producer_finished = false;

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

    producer_finished.store(true);

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
void worker(SafeQ<int> &q, std::atomic<int> &primes, std::atomic<int> &nonprimes, std::atomic<double> &sum,
            std::atomic<int> &consumed_count, std::vector<int> &number_counts)
{
    // implement: use synchronization
    // Note: This part may need some rearranging and rewriting
    // the while loop cannot just check for the size,
    // it has to now wait until the next element can be popped,
    // or it has to terminate if producer has finished and the queue is empty.

    for(;;)
    {
        int num = q.wait_and_pop();
        std::cout << num << std::endl;
        if (producer_finished.load() && q.empty()) {
            break;
        }

        const bool is_prime = kernel(num) == 1;

        consumed_count++;
        number_counts[num % 10]++;

        // sum += num;

        if (is_prime) primes++;
        else nonprimes++;
    }
}

int main(int argc, char **argv)
{
    int num_threads = 10; // you can change this default to thread::hardware_concurrency()
    bool no_exec_times = false, only_exec_times = false;; // reporting of time measurements
    std::string filename = "input.txt";
    parse_args(argc, argv, num_threads, filename, no_exec_times, only_exec_times);

    // The actual code
    std::atomic<int> primes = 0, nonprimes = 0, count = 0;
    std::atomic<int> consumed_count = 0;
    std::atomic<double> mean = 0.0, sum = 0.0;

    // vector for storing numbers ending with different digits (0-9)
    std::vector<int> number_counts(10, 0);

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
                             std::ref(consumed_count), std::ref(number_counts));
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