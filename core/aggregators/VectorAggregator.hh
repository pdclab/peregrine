#ifndef VECTOR_AGGREGATOR_HH
#define VECTOR_AGGREGATOR_HH

#include "../Options.hh"
#include "../Barrier.hh"
#include <vector>

namespace Peregrine
{
  template <typename AggValueT>
  struct VecAggItem
  {
    std::vector<AggValueT> *v;
    bool fresh;
  };
  
  template <typename AggValueT, OnTheFlyOption onthefly, StoppableOption stoppable, typename ViewFunc>
  struct VecAggHandle;
  
  template <typename AggValueT, OnTheFlyOption onthefly, StoppableOption stoppable, typename ViewFunc>
  struct VecAggregator
  {
    using ViewType = decltype(std::declval<ViewFunc>()(std::declval<AggValueT>()));
    using AggHandle = VecAggHandle<AggValueT, onthefly, stoppable, ViewFunc>;
  
    VecAggregator(uint32_t nworkers, ViewFunc &vf)
      : VEC_AGG_OFFSET(Context::data_graph->get_label_range().first),
        VEC_AGG_SIZE(Context::data_graph->get_label_range().second - VEC_AGG_OFFSET),
        values(nworkers),
        handles(nworkers),
        flag({false, false}),
        viewer(vf),
        latest_result(VEC_AGG_SIZE)
    {}
  
    VecAggregator(VecAggregator &) = delete;
  
    const uint32_t VEC_AGG_OFFSET;
    const uint32_t VEC_AGG_SIZE;
    std::vector<AggValueT> global;
    std::vector<std::atomic<VecAggItem<AggValueT>>> values;
    std::vector<AggHandle *> handles;
    std::atomic<flag_t> flag;
    ViewFunc viewer;
    std::vector<std::atomic<ViewType>> latest_result;
  
    bool stale(uint32_t id) const
    {
      return !values[id].load().fresh;
    }
  
    void set_fresh(uint32_t id)
    {
      values[id].store({values[id].load().v, true});
    }
  
    void update()
    {
      flag_t expected = ON();
      if (flag.compare_exchange_strong(expected, WORKING()))
      {
        update_unchecked();
        flag_t e = WORKING();
        flag.compare_exchange_strong(e, ON());
      }
    }
  
    void update_unchecked()
    {
      for (auto &val : values)
      {
        if (val.load().fresh)
        {
          std::vector<AggValueT> *vec = val.load().v;
          for (uint32_t i = 0; i < vec->size(); ++i)
          {
            global[i] += (*vec)[i];
            latest_result[i].store(viewer(global[i]));
          }
          val.store({vec, false});
        }
      }
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
  
    void reset()
    {
      global.clear();
      global.resize(VEC_AGG_SIZE);
      for (auto &item : latest_result)
      {
        item.store(ViewType());
      }
  
      // cmp_exchg isn't really necessary but
      flag_t expected = ON();
      flag.compare_exchange_strong(expected, OFF());
    }
  
    void register_handle(uint32_t id, AggHandle &ah)
    {
      values[id] = {&ah.other, false};
      handles[id] = &ah;
    }
  };
  
  template <typename AggValueT, OnTheFlyOption onthefly, StoppableOption stoppable, typename ViewFunc>
  struct VecAggHandle
  {
    using ViewType = decltype(std::declval<ViewFunc>()(std::declval<AggValueT>()));
    using Aggregator = VecAggregator<AggValueT, onthefly, stoppable, ViewFunc>;

    const uint32_t VEC_AGG_OFFSET;
    const uint32_t VEC_AGG_SIZE;
  
    std::vector<AggValueT> curr;
    std::vector<AggValueT> other;
  
    uint32_t id;
    Aggregator *agg;
    uint32_t new_label_idx;
    Barrier &barrier;

    VecAggHandle(uint32_t tid, Aggregator *a, Barrier &b)
      : VEC_AGG_OFFSET(Context::data_graph->get_label_range().first),
        VEC_AGG_SIZE(Context::data_graph->get_label_range().second - VEC_AGG_OFFSET),
        curr(VEC_AGG_SIZE),
        other(VEC_AGG_SIZE),
        id(tid),
        agg(a),
        new_label_idx(Context::data_graph->new_label),
        barrier(b)
    {}
  
    // Important! Writing AggValueT &v ruins performance.
    // Consider FSM: AggValueT is Domain, but we map vectors, so each call to
    // map constructs a Domain from the vector. Using auto instead means
    // we use the vector, saving the unnecessary construction.
    void map(const std::vector<uint32_t> &k, const auto &v)
    {
      curr[k[new_label_idx] - VEC_AGG_OFFSET] += v;
    }
  
    void reset()
    {
      curr.clear();
      curr.resize(VEC_AGG_SIZE);
  
      other.clear();
      other.resize(VEC_AGG_SIZE);

      new_label_idx = Context::data_graph->new_label;
    }
  
    ViewType read_value(const std::vector<uint32_t> &k)
    {
      return agg->latest_result[k[new_label_idx] - VEC_AGG_OFFSET].load();
    }
  
    void stop()
    {
      barrier.stopAll();
    }
  
    void submit()
    {
      // if other has been read
      if (agg->stale(id))
      {
        // swap curr and other
        std::swap(curr, other);
        curr.clear();
        curr.resize(VEC_AGG_SIZE);
  
        // set freshness
        agg->set_fresh(id);
      }
    }
  };
}

#endif
