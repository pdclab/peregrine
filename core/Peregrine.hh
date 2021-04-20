#ifndef PEREGRINE_HH
#define PEREGRINE_HH

#include <type_traits>
#include <thread>
#include <chrono>
#include <mutex>
#include <condition_variable>
#include <memory>

#include "Options.hh"
#include "Graph.hh"
#include "PatternGenerator.hh"
#include "PatternMatching.hh"

#define CALL_COUNT_LOOP(L, has_anti_vertices)\
{\
  switch (L)\
  {\
    case Graph::LABELLED:\
      lcount += count_loop<Graph::LABELLED, has_anti_vertices>(dg, cands);\
      break;\
    case Graph::UNLABELLED:\
      lcount += count_loop<Graph::UNLABELLED, has_anti_vertices>(dg, cands);\
      break;\
    case Graph::PARTIALLY_LABELLED:\
      lcount += count_loop<Graph::PARTIALLY_LABELLED, has_anti_vertices>(dg, cands);\
      break;\
    case Graph::DISCOVER_LABELS:\
      lcount += count_loop<Graph::DISCOVER_LABELS, has_anti_vertices>(dg, cands);\
      break;\
  }\
}

#define CALL_MATCH_LOOP(L, has_anti_edges, has_anti_vertices)\
{\
  switch (L)\
  {\
    case Graph::LABELLED:\
      match_loop<Graph::LABELLED, has_anti_edges, has_anti_vertices, OnTheFly, Stoppable>(dg, process, cands, ah);\
      break;\
    case Graph::UNLABELLED:\
      match_loop<Graph::UNLABELLED, has_anti_edges, has_anti_vertices, OnTheFly, Stoppable>(dg, process, cands, ah);\
      break;\
    case Graph::PARTIALLY_LABELLED:\
      match_loop<Graph::PARTIALLY_LABELLED, has_anti_edges, has_anti_vertices, OnTheFly, Stoppable>(dg, process, cands, ah);\
      break;\
    case Graph::DISCOVER_LABELS:\
      match_loop<Graph::DISCOVER_LABELS, has_anti_edges, has_anti_vertices, OnTheFly, Stoppable>(dg, process, cands, ah);\
      break;\
  }\
}

namespace Peregrine
{
  // XXX: construct for each application?
  // so that e.g. gcount isn't accessible from a match()
  // or so task_ctr can't be modified at runtime
  namespace Context
  {
    std::shared_ptr<AnalyzedPattern> current_pattern;
    DataGraph *data_graph;
    std::atomic<uint64_t> task_ctr(0);
    std::atomic<uint64_t> gcount(0);
  }

  struct flag_t { bool on, working; };
  static constexpr flag_t ON() { return {true, false}; }
  static constexpr flag_t WORKING() { return {true, true}; }
  static constexpr flag_t OFF() { return {false, false}; }
}


#include "aggregators/SingleValueAggregator.hh"
#include "aggregators/VectorAggregator.hh"
#include "aggregators/Aggregator.hh"
#include "Barrier.hh"

namespace Peregrine
{
  template <Graph::Labelling L,
    bool has_anti_edges,
    bool has_anti_vertices,
    OnTheFlyOption OnTheFly,
    StoppableOption Stoppable,
    typename Func,
    typename HandleType>
  inline void match_loop(DataGraph *dg, const Func &process, std::vector<std::vector<uint32_t>> &cands, HandleType &a)
  {
    uint32_t vgs_count = dg->get_vgs_count();
    uint32_t num_vertices = dg->get_vertex_count();
    uint64_t num_tasks = num_vertices * vgs_count;

    uint64_t task = 0;
    while ((task = Context::task_ctr.fetch_add(1, std::memory_order_relaxed) + 1) <= num_tasks)
    {
      uint32_t v = (task-1) / vgs_count + 1;
      uint32_t vgsi = task % vgs_count;
      Matcher<has_anti_vertices, Stoppable, decltype(process)> m(dg->rbi, dg, vgsi, cands, process);
      m.template map_into<L, has_anti_edges>(v);

      if constexpr (Stoppable == STOPPABLE)
      {
        pthread_testcancel();
      }

      if constexpr (OnTheFly == ON_THE_FLY)
      {
        a.submit();
      }
    }
  }

  template <Graph::Labelling L, bool has_anti_vertices>
  inline uint64_t count_loop(DataGraph *dg, std::vector<std::vector<uint32_t>> &cands)
  {
    uint32_t vgs_count = dg->get_vgs_count();
    uint32_t num_vertices = dg->get_vertex_count();
    uint64_t num_tasks = num_vertices * vgs_count;

    uint64_t lcount = 0;

    uint64_t task = 0;
    while ((task = Context::task_ctr.fetch_add(1, std::memory_order_relaxed) + 1) <= num_tasks)
    {
      uint32_t v = (task-1) / vgs_count + 1;
      uint32_t vgsi = task % vgs_count;
      Counter<has_anti_vertices> m(dg->rbi, dg, vgsi, cands);
      lcount += m.template map_into<L>(v);
    }

    return lcount;
  }

  void count_worker(unsigned tid, DataGraph *dg, Barrier &b)
  {
    (void)tid; // unused

    // an extra pre-allocated cand vector for scratch space, and one for anti-vertex
    std::vector<std::vector<uint32_t>> cands(dg->rbi.query_graph.num_vertices() + 2);

    while (b.hit())
    {
      Graph::Labelling L = dg->rbi.labelling_type();
      bool has_anti_edges = dg->rbi.has_anti_edges();
      bool has_anti_vertices = !dg->rbi.anti_vertices.empty();

      cands.resize(dg->rbi.query_graph.num_vertices()+2);
      for (auto &cand : cands) {
        cand.clear();
        cand.reserve(10000);
      }

      uint64_t lcount = 0;

      if (has_anti_edges)
      {
        constexpr StoppableOption Stoppable = UNSTOPPABLE;
        constexpr OnTheFlyOption OnTheFly = AT_THE_END;

        const auto process = [&lcount](const CompleteMatch &) { lcount += 1; };
        // dummy
        struct {void submit() {}} ah;

        // TODO anti-edges ruin a lot of optimizations, but not all,
        // is there no way to handle them in Counter?
        if (has_anti_vertices)
        {
          CALL_MATCH_LOOP(L, true, true);
        }
        else
        {
          CALL_MATCH_LOOP(L, true, false);
        }
      }
      else
      {
        if (has_anti_vertices)
        {
          CALL_COUNT_LOOP(L, true);
        }
        else
        {
          CALL_COUNT_LOOP(L, false);
        }
      }

      Context::gcount += lcount;
    }
  }

  template <
    typename AggKeyT,
    typename AggValueT,
    OnTheFlyOption OnTheFly,
    StoppableOption Stoppable,
    typename AggregatorType,
    typename F
  >
  void map_worker(unsigned tid, DataGraph *dg, Barrier &b, AggregatorType &a, F &&p)
  {
    // an extra pre-allocated cand vector for scratch space, and one for anti-vertex
    std::vector<std::vector<uint32_t>> cands(dg->rbi.query_graph.num_vertices() + 2);

    using ViewFunc = decltype(a.viewer);
    MapAggHandle<AggKeyT, AggValueT, OnTheFly, Stoppable, ViewFunc> ah(tid, &a, b);
    a.register_handle(tid, ah);

    while (b.hit())
    {
      ah.reset();
      const auto process = [&ah, &p](const CompleteMatch &cm) { p(ah, cm); };
      Graph::Labelling L = dg->rbi.labelling_type();
      bool has_anti_edges = dg->rbi.has_anti_edges();
      bool has_anti_vertices = !dg->rbi.anti_vertices.empty();

      cands.resize(dg->rbi.query_graph.num_vertices()+2);
      for (auto &cand : cands) {
        cand.clear();
        cand.reserve(10000);
      }

      if (has_anti_edges)
      {
        if (has_anti_vertices)
        {
          CALL_MATCH_LOOP(L, true, true);
        }
        else
        {
          CALL_MATCH_LOOP(L, true, false);
        }
      }
      else
      {
        if (has_anti_vertices)
        {
          CALL_MATCH_LOOP(L, false, true);
        }
        else
        {
          CALL_MATCH_LOOP(L, false, false);
        }
      }
    }
  }


  template <
    typename AggValueT,
    OnTheFlyOption OnTheFly,
    StoppableOption Stoppable,
    typename AggregatorType,
    typename F
  >
  void single_worker(unsigned tid, DataGraph *dg, Barrier &b, AggregatorType &a, F &&p)
  {
    // an extra pre-allocated cand vector for scratch space, and one for anti-vertex
    std::vector<std::vector<uint32_t>> cands(dg->rbi.query_graph.num_vertices() + 2);

    using ViewFunc = decltype(a.viewer);
    SVAggHandle<AggValueT, OnTheFly, Stoppable, ViewFunc> ah(tid, &a, b);
    a.register_handle(tid, ah);

    while (b.hit())
    {
      ah.reset();
      const auto process = [&ah, &p](const CompleteMatch &cm) { p(ah, cm); };

      Graph::Labelling L = dg->rbi.labelling_type();
      bool has_anti_edges = dg->rbi.has_anti_edges();
      bool has_anti_vertices = !dg->rbi.anti_vertices.empty();

      cands.resize(dg->rbi.query_graph.num_vertices()+2);
      for (auto &cand : cands) {
        cand.clear();
        cand.reserve(10000);
      }

      if (has_anti_edges)
      {
        if (has_anti_vertices)
        {
          CALL_MATCH_LOOP(L, true, true);
        }
        else
        {
          CALL_MATCH_LOOP(L, true, false);
        }
      }
      else
      {
        if (has_anti_vertices)
        {
          CALL_MATCH_LOOP(L, false, true);
        }
        else
        {
          CALL_MATCH_LOOP(L, false, false);
        }
      }
    }
  }

  template <
    typename AggValueT,
    OnTheFlyOption OnTheFly,
    StoppableOption Stoppable,
    typename AggregatorType,
    typename F
  >
  void vector_worker(unsigned tid, DataGraph *dg, Barrier &b, AggregatorType &a, F &&p)
  {
    // an extra pre-allocated cand vector for scratch space, and one for anti-vertex
    std::vector<std::vector<uint32_t>> cands(dg->rbi.query_graph.num_vertices() + 2);

    using ViewFunc = decltype(a.viewer);
    VecAggHandle<AggValueT, OnTheFly, Stoppable, ViewFunc> ah(tid, &a, b);
    a.register_handle(tid, ah);

    while (b.hit())
    {
      ah.reset();
      const auto process = [&ah, &p](const CompleteMatch &cm) { p(ah, cm); };

      Graph::Labelling L = dg->rbi.labelling_type();
      bool has_anti_edges = dg->rbi.has_anti_edges();
      bool has_anti_vertices = !dg->rbi.anti_vertices.empty();

      cands.resize(dg->rbi.query_graph.num_vertices()+2);
      for (auto &cand : cands) {
        cand.clear();
        cand.reserve(10000);
      }

      if (has_anti_edges)
      {
        if (has_anti_vertices)
        {
          CALL_MATCH_LOOP(L, true, true);
        }
        else
        {
          CALL_MATCH_LOOP(L, true, false);
        }
      }
      else
      {
        if (has_anti_vertices)
        {
          CALL_MATCH_LOOP(L, false, true);
        }
        else
        {
          CALL_MATCH_LOOP(L, false, false);
        }
      }
    }
  }

  template <typename AggregatorType>
  void aggregator_thread(Barrier &barrier, AggregatorType &agg)
  {
    using namespace std::chrono_literals;

    while (!barrier.finished())
    {
      agg.update();
      std::this_thread::sleep_for(300ms);
    }
  }

  template <typename T>
  struct trivial_wrapper
  {
    trivial_wrapper() : val() {}
    trivial_wrapper(T v) : val(v) {}
    trivial_wrapper<T> &operator+=(const trivial_wrapper<T> &other) { val += other.val; return *this; }
    void reset() { val = T(); }
    T val;
  };

  template <typename T>
  T default_viewer(T &&v) { return v; }

  template <
    typename AggKeyT,
    typename GivenAggValueT,
    OnTheFlyOption OnTheFly,
    StoppableOption Stoppable,
    typename DataGraphT,
    typename PF,
    typename VF = decltype(default_viewer<GivenAggValueT>)
  >
  std::vector<std::pair<SmallGraph, decltype(std::declval<VF>()(std::declval<GivenAggValueT>()))>>
  match(DataGraphT &&data_graph,
      const std::vector<SmallGraph> &patterns,
      size_t nworkers,
      PF &&process,
      VF viewer = default_viewer<GivenAggValueT>)
  {
    if (patterns.empty())
    {
      return {};
    }

    // automatically wrap trivial types so they have .reset() etc
    constexpr bool should_be_wrapped = std::is_trivial<GivenAggValueT>::value;
    using AggValueT = typename std::conditional<should_be_wrapped,
      trivial_wrapper<GivenAggValueT>, GivenAggValueT>::type;
    auto view = [&viewer](auto &&v)
    {
      if constexpr (should_be_wrapped)
      {
        return viewer(std::move(v.val));
      }
      else
      {
        return viewer(v);
      }
    };

    if constexpr (std::is_same_v<std::decay_t<DataGraphT>, DataGraph>)
    {
      Context::data_graph = &data_graph;
    }
    else if constexpr (std::is_same_v<std::decay_t<DataGraphT>, DataGraph *>)
    {
      Context::data_graph = data_graph;
    }
    else
    {
      Context::data_graph = new DataGraph(data_graph);
      utils::Log{} << "Finished reading datagraph: |V| = " << Context::data_graph->get_vertex_count()
                << " |E| = " << Context::data_graph->get_edge_count()
                << "\n";
    }

    std::vector<std::pair<SmallGraph, decltype(std::declval<VF>()(std::declval<GivenAggValueT>()))>> result;

    // optimize AggKeyT == Pattern
    if constexpr (std::is_same_v<AggKeyT, Pattern>)
    {
      Context::data_graph->set_known_labels(patterns);

      std::vector<SmallGraph> single;
      std::vector<SmallGraph> vector;
      std::vector<SmallGraph> multi;

      for (const auto &p : patterns)
      {
        Graph::Labelling l = p.get_labelling();
        switch (l)
        {
          case Graph::LABELLED:
          case Graph::UNLABELLED:
            single.emplace_back(p);
            break;
          case Graph::PARTIALLY_LABELLED:
            vector.emplace_back(p);
            break;
          case Graph::DISCOVER_LABELS:
            multi.emplace_back(p);
            break;
        }
      }

      if (!single.empty()
          && std::is_integral_v<GivenAggValueT>
          && Stoppable == UNSTOPPABLE && OnTheFly == AT_THE_END)
      {
        utils::Log{}
          << "WARN: If you are counting, Peregrine::count() is much faster!"
          << "\n";
      }

      result = match_single<AggValueT, OnTheFly, Stoppable>(process, view, nworkers, single);
      auto vector_result = match_vector<AggValueT, OnTheFly, Stoppable>(process, view, nworkers, vector);
      auto multi_result = match_multi<AggKeyT, AggValueT, OnTheFly, Stoppable>(process, view, nworkers, multi);

      result.insert(result.end(), vector_result.begin(), vector_result.end());
      result.insert(result.end(), multi_result.begin(), multi_result.end());

      return result;
    }
    else
    {
      result = match_multi<AggKeyT, AggValueT, OnTheFly, Stoppable>(process, view, nworkers, patterns);
    }

    if constexpr (!std::is_same_v<std::decay_t<DataGraphT>, DataGraph> && !std::is_same_v<std::decay_t<DataGraphT>, DataGraph *>)
    {
      delete Context::data_graph;
    }

    return result;
  }

  template <typename AggKeyT, typename AggValueT, OnTheFlyOption OnTheFly, StoppableOption Stoppable, typename PF, typename VF>
  std::vector<std::pair<SmallGraph, decltype(std::declval<VF>()(std::declval<AggValueT>()))>>
  match_multi
  (PF &&process, VF &&viewer, size_t nworkers, const std::vector<SmallGraph> &patterns)
  {
    std::vector<std::pair<SmallGraph, decltype(std::declval<VF>()(std::declval<AggValueT>()))>> results;

    if (patterns.empty()) return results;

    // initialize
    Barrier barrier(nworkers);
    std::vector<std::thread> pool;
    DataGraph *dg(Context::data_graph);
    dg->set_rbi(patterns.front());

    Context::current_pattern = std::make_shared<AnalyzedPattern>(AnalyzedPattern(dg->rbi));

    MapAggregator<AggKeyT, AggValueT, OnTheFly, Stoppable, decltype(viewer)> aggregator(nworkers, viewer);

    for (uint8_t i = 0; i < nworkers; ++i)
    {
      //auto &ah = aggregator.get_handle(i);
      pool.emplace_back(map_worker<
            AggKeyT,
            AggValueT,
            OnTheFly,
            Stoppable,
            decltype(aggregator),
            PF
          >,
          i,
          dg,
          std::ref(barrier),
          std::ref(aggregator),
          std::ref(process));
    }

    std::thread agg_thread;
    if constexpr (OnTheFly == ON_THE_FLY)
    {
      agg_thread = std::thread(aggregator_thread<decltype(aggregator)>, std::ref(barrier), std::ref(aggregator));
    }

    // make sure the threads are all running
    barrier.join();

    auto t1 = utils::get_timestamp();
    for (const auto &p : patterns)
    {
      // reset state
      Context::task_ctr = 0;

      // set new pattern
      dg->set_rbi(p);
      Context::current_pattern = std::make_shared<AnalyzedPattern>(AnalyzedPattern(dg->rbi));
      // prepare handles for the next pattern
      aggregator.reset();

      // begin matching
      barrier.release();

      // sleep until matching finished
      bool called_stop = barrier.join();

      aggregator.get_result();

      if constexpr (Stoppable == STOPPABLE)
      {
        if (called_stop)
        {
          // cancel
          for (auto &th : pool)
          {
            pthread_cancel(th.native_handle());
          }

          // wait for them all to end
          for (auto &th : pool)
          {
            th.join();
          }

          pool.clear();
          barrier.reset();

          // restart workers
          for (uint8_t i = 0; i < nworkers; ++i)
          {
            pool.emplace_back(map_worker<
                  AggKeyT,
                  AggValueT,
                  OnTheFly,
                  Stoppable,
                  decltype(aggregator),
                  PF
                >,
                i,
                dg,
                std::ref(barrier),
                std::ref(aggregator),
                std::ref(process));
          }

          barrier.join();
        }
      }

      for (auto &[k, v] : aggregator.latest_result)
      {
        results.emplace_back(SmallGraph(p, k), v);
      }
    }
    auto t2 = utils::get_timestamp();

    barrier.finish();
    for (auto &th : pool)
    {
      th.join();
    }

    if constexpr (OnTheFly == ON_THE_FLY)
    {
      agg_thread.join();
    }

    utils::Log{} << "-------" << "\n";
    utils::Log{} << "all patterns finished after " << (t2-t1)/1e6 << "s" << "\n";

    return results;
  }

  template <typename AggValueT, OnTheFlyOption OnTheFly, StoppableOption Stoppable, typename PF, typename VF>
  std::vector<std::pair<SmallGraph, decltype(std::declval<VF>()(std::declval<AggValueT>()))>>
  match_single
  (PF &&process, VF &&viewer, size_t nworkers, const std::vector<SmallGraph> &patterns)
  {
    std::vector<std::pair<SmallGraph, decltype(std::declval<VF>()(std::declval<AggValueT>()))>> results;

    if (patterns.empty()) return results;

    // initialize
    Barrier barrier(nworkers);
    std::vector<std::thread> pool;
    DataGraph *dg(Context::data_graph);
    dg->set_rbi(patterns.front());

    Context::current_pattern = std::make_shared<AnalyzedPattern>(AnalyzedPattern(dg->rbi));

    SVAggregator<AggValueT, OnTheFly, Stoppable, decltype(viewer)> aggregator(nworkers, viewer);

    for (uint8_t i = 0; i < nworkers; ++i)
    {
      pool.emplace_back(single_worker<
            AggValueT,
            OnTheFly,
            Stoppable,
            decltype(aggregator),
            PF
          >,
          i,
          dg,
          std::ref(barrier),
          std::ref(aggregator),
          std::ref(process));
    }

    std::thread agg_thread;
    if constexpr (OnTheFly == ON_THE_FLY)
    {
      agg_thread = std::thread(aggregator_thread<decltype(aggregator)>, std::ref(barrier), std::ref(aggregator));
    }

    // make sure the threads are all running
    barrier.join();

    auto t1 = utils::get_timestamp();
    for (const auto &p : patterns)
    {
      // reset state
      Context::task_ctr = 0;

      // set new pattern
      dg->set_rbi(p);
      Context::current_pattern = std::make_shared<AnalyzedPattern>(AnalyzedPattern(dg->rbi));
      // prepare handles for the next pattern
      aggregator.reset();

      // begin matching
      barrier.release();

      // sleep until matching finished
      bool called_stop = barrier.join();

      // need to get thread-local values before killing threads
      aggregator.get_result();
      results.emplace_back(p, aggregator.latest_result.load());

      if constexpr (Stoppable == STOPPABLE)
      {
        if (called_stop)
        {
          // cancel
          for (auto &th : pool)
          {
            pthread_cancel(th.native_handle());
          }

          // wait for them all to end
          for (auto &th : pool)
          {
            th.join();
          }

          pool.clear();
          barrier.reset();

          // restart workers
          for (uint8_t i = 0; i < nworkers; ++i)
          {
            pool.emplace_back(single_worker<
                  AggValueT,
                  OnTheFly,
                  Stoppable,
                  decltype(aggregator),
                  PF
                >,
                i,
                dg,
                std::ref(barrier),
                std::ref(aggregator),
                std::ref(process));
          }

          barrier.join();
        }
      }
    }
    auto t2 = utils::get_timestamp();

    barrier.finish();
    for (auto &th : pool)
    {
      th.join();
    }

    if constexpr (OnTheFly == ON_THE_FLY)
    {
      agg_thread.join();
    }

    utils::Log{} << "-------" << "\n";
    utils::Log{} << "all patterns finished after " << (t2-t1)/1e6 << "s" << "\n";

    return results;
  }

  template <typename AggValueT, OnTheFlyOption OnTheFly, StoppableOption Stoppable, typename PF, typename VF>
  std::vector<std::pair<SmallGraph, decltype(std::declval<VF>()(std::declval<AggValueT>()))>>
  match_vector
  (PF &&process, VF &&viewer, size_t nworkers, const std::vector<SmallGraph> &patterns)
  {
    std::vector<std::pair<SmallGraph, decltype(std::declval<VF>()(std::declval<AggValueT>()))>> results;

    if (patterns.empty()) return results;

    // initialize
    Barrier barrier(nworkers);
    std::vector<std::thread> pool;
    DataGraph *dg(Context::data_graph);
    dg->set_rbi(patterns.front());

    Context::current_pattern = std::make_shared<AnalyzedPattern>(AnalyzedPattern(dg->rbi));

    VecAggregator<AggValueT, OnTheFly, Stoppable, decltype(viewer)> aggregator(nworkers, viewer);

    for (uint8_t i = 0; i < nworkers; ++i)
    {
      //auto &ah = aggregator.get_handle(i);
      pool.emplace_back(vector_worker<
            AggValueT,
            OnTheFly,
            Stoppable,
            decltype(aggregator),
            PF
          >,
          i,
          dg,
          std::ref(barrier),
          std::ref(aggregator),
          std::ref(process));
    }

    std::thread agg_thread;
    if constexpr (OnTheFly == ON_THE_FLY)
    {
      agg_thread = std::thread(aggregator_thread<decltype(aggregator)>, std::ref(barrier), std::ref(aggregator));
    }

    // make sure the threads are all running
    barrier.join();

    auto t1 = utils::get_timestamp();
    for (const auto &p : patterns)
    {
      // reset state
      Context::task_ctr = 0;

      // set new pattern
      dg->set_rbi(p);
      Context::current_pattern = std::make_shared<AnalyzedPattern>(AnalyzedPattern(dg->rbi));
      // prepare handles for the next pattern
      aggregator.reset();

      // begin matching
      barrier.release();

      // sleep until matching finished
      bool called_stop = barrier.join();

      aggregator.get_result();

      if constexpr (Stoppable == STOPPABLE)
      {
        if (called_stop)
        {
          // cancel
          for (auto &th : pool)
          {
            pthread_cancel(th.native_handle());
          }

          // wait for them all to end
          for (auto &th : pool)
          {
            th.join();
          }

          pool.clear();
          barrier.reset();

          // restart workers
          for (uint8_t i = 0; i < nworkers; ++i)
          {
            pool.emplace_back(vector_worker<
                  AggValueT,
                  OnTheFly,
                  Stoppable,
                  decltype(aggregator),
                  PF
                >,
                i,
                dg,
                std::ref(barrier),
                std::ref(aggregator),
                std::ref(process));
          }

          barrier.join();
        }
      }

      std::vector<uint32_t> ls(p.get_labels().cbegin(), p.get_labels().cend());
      uint32_t pl = dg->new_label;
      uint32_t l = 0;
      for (auto &m : aggregator.latest_result)
      {
        ls[pl] = aggregator.VEC_AGG_OFFSET + l;
        results.emplace_back(SmallGraph(p, ls), m.load());
        l += 1;
      }
    }
    auto t2 = utils::get_timestamp();

    barrier.finish();
    for (auto &th : pool)
    {
      th.join();
    }

    if constexpr (OnTheFly == ON_THE_FLY)
    {
      agg_thread.join();
    }

    utils::Log{} << "-------" << "\n";
    utils::Log{} << "all patterns finished after " << (t2-t1)/1e6 << "s" << "\n";

    return results;
  }

  // for each pattern, calculate the vertex-based count
  std::vector<std::pair<SmallGraph, uint64_t>> convert_counts(std::vector<std::pair<SmallGraph, uint64_t>> edge_based, const std::vector<SmallGraph> &original_patterns)
  {
    std::vector<std::pair<SmallGraph, uint64_t>> vbased(edge_based.size());

    for (int32_t i = edge_based.size()-1; i >= 0; --i) {
      uint64_t count = edge_based[i].second;
      for (uint32_t j = i+1; j < edge_based.size(); ++j) {
        // mapping edge_based[i].first into edge_based[j].first
        uint32_t n = num_mappings(edge_based[j].first, edge_based[i].first);
        uint64_t inc = n * vbased[j].second;
        count -= inc;
      }
      vbased[i] = {original_patterns[i], count};
    }

    return vbased;
  }

  template <typename DataGraphT>
  std::vector<std::pair<SmallGraph, uint64_t>>
  count(DataGraphT &&data_graph, const std::vector<SmallGraph> &patterns, size_t nworkers)
  {
    // initialize
    std::vector<std::pair<SmallGraph, uint64_t>> results;
    if (patterns.empty()) return results;

    // optimize if all unlabelled vertex-induced patterns of a certain size
    // TODO: if a subset is all unlabelled vertex-induced patterns of a certain
    // size it can be optimized too
    uint32_t sz = patterns.front().num_vertices();
    auto is_same_size = [&sz](const SmallGraph &p) {
        return p.num_vertices() == sz && p.num_anti_vertices() == 0;
      };
    auto is_unlabelled = [&sz](const SmallGraph &p) {
        return p.get_labelling() == Graph::UNLABELLED;
      };
    auto is_vinduced = [](const SmallGraph &p) {
        uint32_t m = p.num_anti_edges() + p.num_true_edges();
        uint32_t n = p.num_vertices();
        return m == (n*(n-1))/2;
      };
    uint32_t num_possible_topologies[] = {
      0,
      1,
      1,
      2, // size 3
      6, // size 4
      21, // size 5
      112, // size 6
      853, // size 7
      11117, // size 8
      261080, // size 9
    };

    bool must_convert_counts = false;
    std::vector<SmallGraph> new_patterns;
    if (std::all_of(patterns.cbegin(), patterns.cend(), is_same_size)
        && std::all_of(patterns.cbegin(), patterns.cend(), is_unlabelled)
        && std::all_of(patterns.cbegin(), patterns.cend(), is_vinduced)
        && (sz < 10 && patterns.size() == num_possible_topologies[sz]))
    {
      must_convert_counts = true;
      new_patterns = PatternGenerator::all(sz, PatternGenerator::VERTEX_BASED, PatternGenerator::EXCLUDE_ANTI_EDGES);
    }
    else
    {
      new_patterns.assign(patterns.cbegin(), patterns.cend());
    }

    Barrier barrier(nworkers);
    std::vector<std::thread> pool;

    DataGraph *dg;
    if constexpr (std::is_same_v<std::decay_t<DataGraphT>, DataGraph>)
    {
      dg = &data_graph;
    }
    else if constexpr (std::is_same_v<std::decay_t<DataGraphT>, DataGraph *>)
    {
      dg = data_graph;
    }
    else
    {
      dg = new DataGraph(data_graph);
      utils::Log{} << "Finished reading datagraph: |V| = " << dg->get_vertex_count()
                << " |E| = " << dg->get_edge_count()
                << "\n";
    }

    dg->set_rbi(new_patterns.front());
    dg->set_known_labels(new_patterns);

    for (uint8_t i = 0; i < nworkers; ++i)
    {
      pool.emplace_back(count_worker,
          i,
          dg,
          std::ref(barrier));
    }

    // make sure the threads are all running
    barrier.join();

    auto t1 = utils::get_timestamp();
    for (const auto &p : new_patterns)
    {
      // reset state
      Context::task_ctr = 0;
      Context::gcount = 0;

      // set new pattern
      dg->set_rbi(p);

      // begin matching
      barrier.release();

      // sleep until matching finished
      barrier.join();

      // get counts
      uint64_t global_count = Context::gcount;
      results.emplace_back(p, global_count);
    }
    auto t2 = utils::get_timestamp();

    barrier.finish();
    for (auto &th : pool)
    {
      th.join();
    }

    if (must_convert_counts)
    {
      results = convert_counts(results, patterns);
    }

    if constexpr (!std::is_same_v<std::decay_t<DataGraphT>, DataGraph> && !std::is_same_v<std::decay_t<DataGraphT>, DataGraph *>)
    {
      delete dg;
    }

    utils::Log{} << "-------" << "\n";
    utils::Log{} << "all patterns finished after " << (t2-t1)/1e6 << "s" << "\n";


    return results;
  }
} // namespace Peregrine

#endif
