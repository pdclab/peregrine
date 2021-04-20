#ifndef BARRIER_HH
#define BARRIER_HH

namespace Peregrine
{
  enum WORKER_STATE { WAIT, GO, DONE };
  class Barrier
  {
    public:
      Barrier(size_t n) : s(WAIT), n_workers(n), workers_waiting(0)
      {
      }

      bool hit()
      {
        // wait for other threads to catch the GO signal
        while (s.load(std::memory_order_relaxed) == GO);
        // at this point, s == WAIT

        // indicate you've hit the barrier
        uint32_t n = workers_waiting.fetch_add(1, std::memory_order_relaxed);
        if (n+1 == n_workers)
        {
          // notify coordinator
          {
            std::lock_guard l(mtx);
            all_workers_waiting = true;
          }
          workers_cv.notify_all();
        }

        // wait until coordinator says what to do
        while (s.load(std::memory_order_relaxed) == WAIT);

        if (s.load(std::memory_order_relaxed) == GO)
        {
          // indicate you've left the barrier
          size_t waiting = workers_waiting.fetch_sub(1, std::memory_order_relaxed) - 1;
          if (waiting == 0)
          {
            // signal that all workers have begun working
            s.store(WAIT, std::memory_order_relaxed);
          }

          return true;
        }
        else // s == DONE
        {
          return false;
        }
      }

      bool join()
      {
        std::unique_lock l(mtx);
        workers_cv.wait(l, [this]() { return all_workers_waiting; });
        all_workers_waiting = false;
        return kill_all;
      }

      void stopAll()
      {
        {
          std::unique_lock l(mtx);
          kill_all = true;
          all_workers_waiting = true;
        }
        workers_cv.notify_all();
        pthread_exit(0);
      }

      void reset()
      {
        s = WAIT;
        kill_all = false;
        all_workers_waiting = false;
        workers_waiting = 0;
      }

      // releases all threads
      void release()
      {
        s = GO;
      }

      // instructs all threads to exit
      void finish()
      {
        s = DONE;
      }

      bool finished()
      {
        return s.load(std::memory_order_relaxed) == DONE;
      }

    private:
      std::atomic<WORKER_STATE> s;
      const size_t n_workers;
      bool all_workers_waiting = false;
      bool kill_all = false;
      std::atomic<size_t> workers_waiting;
      std::mutex mtx;
      std::condition_variable workers_cv;
  };
}

#endif
