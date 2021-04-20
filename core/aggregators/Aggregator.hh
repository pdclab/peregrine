#ifndef AGGREGATOR_HH
#define AGGREGATOR_HH

#include "../Options.hh"
#include "../Barrier.hh"

#include "tbb/concurrent_unordered_map.h"
#include <unordered_map>

namespace std
{
  template <>
  struct hash<Peregrine::Pattern>
  {
    std::size_t operator()(const Peregrine::Pattern &k) const
    {
      std::size_t seed = k.size();

      for (auto &i : k) {
        seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };
}

namespace Peregrine
{
  template <typename AggKeyT, typename AggValueT>
  struct MapAggItem
  {
    std::unordered_map<AggKeyT, AggValueT> *v;
    bool fresh;
  };
  
  template <typename AggKeyT, typename AggValueT, OnTheFlyOption onthefly, StoppableOption stoppable, typename ViewFunc>
  class MapAggHandle;
  
  template <typename AggKeyT, typename AggValueT, OnTheFlyOption onthefly, StoppableOption stoppable, typename ViewFunc>
  struct MapAggregator
  {
    using ViewType = decltype(std::declval<ViewFunc>()(std::declval<AggValueT>()));
    using AggHandle = MapAggHandle<AggKeyT, AggValueT, onthefly, stoppable, ViewFunc>;
  
    MapAggregator(uint32_t nworkers, ViewFunc &vf)
      : values(nworkers),
        handles(nworkers),
        flag({false, false}),
        viewer(vf)
    {}
  
    MapAggregator(MapAggregator &) = delete;
    ~MapAggregator() { for (auto handle : handles) delete handle; }
  
    std::unordered_map<AggKeyT, AggValueT> global;
    std::vector<std::atomic<MapAggItem<AggKeyT, AggValueT>>> values;
    std::vector<AggHandle *> handles;
    std::atomic<flag_t> flag;
    ViewFunc viewer;
    tbb::concurrent_unordered_map<AggKeyT, ViewType, std::hash<AggKeyT>> latest_result;
  
    bool stale(uint32_t id) const
    {
      return !values[id].load().fresh;
    }
  
    void set_fresh(uint32_t id)
    {
      values[id].store({values[id].load().v, true});
    }
  
    void get_result()
    {
      // wait for aggregator thread to stop
      flag_t expected = ON();
      while (!flag.compare_exchange_weak(expected, OFF()));
  
      update_unchecked();
      for (auto &handle : handles)
      {
        handle->submit();
      }
      update_unchecked();
    }
  
    void update()
    {
      flag_t expected = ON();
      if (flag.compare_exchange_strong(expected, WORKING()))
      {
        update_unchecked();
        flag_t e = WORKING();
        assert(flag.compare_exchange_strong(e, ON()));
      }
    }
  
    void update_unchecked()
    {
      for (auto &val : values)
      {
        if (val.load().fresh)
        {
          std::unordered_map<AggKeyT, AggValueT> *map = val.load().v;
          //for (uint32_t i = 0; i < vec->size(); ++i)
          for (const auto &[i, v] : *map)
          {
            global[i] += v;
            latest_result[i] = viewer(global[i]);
          }
          val.store({map, false});
        }
      }
    }
  
    void reset()
    {
      global.clear();
      latest_result.clear();
  
      // cmp_exchg isn't really necessary but
      flag_t expected = OFF();
      flag.compare_exchange_strong(expected, ON());
    }
  
    void register_handle(uint32_t id, AggHandle *ah)
    {
      values[id] = {&ah->other, false};
      handles[id] = ah;
    }
  };
  
  template <typename AggKeyT, typename AggValueT, OnTheFlyOption onthefly, StoppableOption stoppable, typename ViewFunc>
  class MapAggHandle
  {
    private:
      using ViewType = decltype(std::declval<ViewFunc>()(std::declval<AggValueT>()));
      using Aggregator = MapAggregator<AggKeyT, AggValueT, onthefly, stoppable, ViewFunc>;
      friend Aggregator;
  
      std::unordered_map<AggKeyT, AggValueT> curr;
      std::unordered_map<AggKeyT, AggValueT> other;
  
      uint32_t id;
      Aggregator *agg;
      Barrier &barrier;
  
    public:
      MapAggHandle(uint32_t tid, Aggregator *a, Barrier &b)
        : id(tid),
          agg(a),
          barrier(b)
      {}
  
      // Important! Writing AggValueT &v ruins performance.
      // Consider FSM: AggValueT is Domain, but we map vectors, so each call to
      // map constructs a Domain from the vector. Using auto instead means
      // we use the vector, saving the unnecessary construction.
      void map(const auto &k, const auto &v)
      {
        curr[k] += v;
      }
  
      void reset()
      {
        curr.clear();
        other.clear();
      }
  
      ViewType read_value(const auto &k)
      {
        return agg->latest_result[k];
      }
  
      void stop()
      {
        if constexpr (stoppable == STOPPABLE)
        {
          barrier.stopAll();
        }
        else
        {
          // XXX: I would like to get a nice error message like this at compile time
          throw std::runtime_error("stop() called on handle to unstoppable aggregator.\n\
              If you want early termination, make sure you call match with the STOPPABLE template argument.\n\
              E.g. `match<Pattern, long, ON_THE_FLY, STOPPABLE>(...)'");
        }
      }
  
      void submit()
      {
        // if other has been read
        if (agg && agg->stale(id))
        {
          // swap curr and other
          std::swap(curr, other);
          curr.clear();
  
          // set freshness
          agg->set_fresh(id);
        }
      }
  };
}

#endif
