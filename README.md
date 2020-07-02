<img src="https://user-images.githubusercontent.com/47194770/78506209-e73ca880-772c-11ea-9fc1-790713527bc9.png" alt="Peregrine" width="70%">

# Peregrine: A Pattern-Aware Graph Mining System

![tests](https://github.com/pdclab/peregrine/workflows/tests/badge.svg)

Peregrine is an efficient, single-machine system for performing data mining tasks on large graphs. Some graph mining applications include:
* Finding frequent subgraphs
* Generating the motif/graphlet distribution
* Finding all occurrences of a subgraph

Peregrine is highly programmable, so you can easily develop your own graph mining applications using its novel, declarative, graph-pattern-centric API.
To write a Peregrine program, you describe which graph patterns you are interested in mining, and what to do with each occurrence of those patterns. You provide the _what_ and the runtime handles the _how_.

For full details, you can read our paper published in [EuroSys 2020](https://dl.acm.org/doi/abs/10.1145/3342195.3387548) or the longer version on [arXiv](https://arxiv.org/abs/2004.02369).

**TL;DR:** compared to other state-of-the-art open-source graph mining systems, Peregrine:
* executes up to 700x faster
* consumes up to 100x less memory
* scales to 100x larger data sets
* on 8x fewer machines
* with a simpler, more expressive API

## Table of Contents
1. [Quick start](#1-quick-start)
2. [Writing your own programs](#2-writing-your-own-programs)
3. [Data graphs](#3-data-graphs)
4. [Reproducing our EuroSys 2020 paper results](#4-reproducing-our-eurosys-2020-paper-results)
5. [Contributing](#5-contributing)
6. [Acknowledgements](#6-acknowledgements)
7. [Resources](#7-resources)

## 1. Quick start

Peregrine has been tested on Ubuntu 18.04 and Arch Linux but should work on any
POSIX-y OS. It requires C++20 support (GCC version >= 9.2.1). Additionally, the
tests require [UnitTest++](https://github.com/unittest-cpp/unittest-cpp).

To build Peregrine:

```
$ git clone https://github.com/pdclab/Peregrine.git
$ cd Peregrine
$ source tbb2020/bin/tbbvars.sh intel64
$ make -j
$ bin/test
```

Several sample applications, query patterns, and a sample dataset are released with the code. Calling any of the applications without arguments will show you what they expect:

```
$ bin/count
USAGE: bin/count <data graph> <pattern | #-motifs | #-clique> [# threads]
```

These applications print their results in `<pattern>: <aggregation value>` format.

For example, motif-counting:

```
$ bin/count data/citeseer 3-motifs 8
Counting 3-motifs
Finished reading datagraph: |V| = 3264 |E| = 4536
[...]
All patterns finished after 0.030265s
[2-3][1-3]: 23380
[1-2][1-3][2-3]: 1166
```

The string `[2-3][1-3]` encodes the pattern consisting of edges `(1, 3), (2, 3)`, and 23380 is the number of unique occurrences Peregrine found in the citeseer graph.

Other applications give similar output:

```
$ bin/count data/citeseer 4-clique 8
[...]
All patterns finished after 0.030265s
[3-4][1-2][1-3][1-4][2-3][2-4]: 255
$
$ bin/count data/citeseer query/p1.graph 8
[...]
All patterns finished after 0.003368s
[3-4][1-2][1-3][1-4][2-3]: 3730
```

FSM provides support values instead of counts:

```
$ bin/fsm data/citeseer 3 300 8 # size 3 FSM with support 300
[...]
Frequent patterns:
[1,0-2,0][1,0-3,0][2,0-4,0]: 303
[1,1-2,1][1,1-3,1][2,1-4,1]: 335
Finished in 0.078629s
```

The existence-query application simply states whether the desired pattern exists or not:

```
$ bin/existence-query data/citeseer 14-clique 8
[...]
All patterns finished after 0.005509s
[pattern omitted due to length] doesn't exist in data/citeseer
```

## 2. Writing your own programs

In Peregrine's programming model, you provide a data graph, a set of patterns
you're interested in, and a callback the system will apply to each occurrence
of these patterns in your data graph. We present a brief overview of the API
here, beginning with constructing patterns.

For all of the following code snippets, assume we are `using namespace Peregrine`.

We have not released support for directed graphs yet; the code currently
assumes all graphs are undirected.

### 2.1 Constructing patterns directly

Pattern graphs can constructed in two ways using the `SmallGraph` data structure:

#### 2.1.1 Statically, from a file

Given a file in the following edge-list format:

```
<vertex-id> [label] <vertex-id> [label]
```

where the `label`'s are optional 32-bit integers. To indicate a vertex is
unlabelled in a partially-labelled pattern, assign it label `-1`. To indicate
an anti-edge, any extra integer can be placed at the end of
the line.

For example, a triangle:

```
1 2
1 3
2 3
```

A vertex-induced 3-star, notice the last edge is an anti-edge:

```
1 2
1 3
2 3 1
```

A partially-labelled triangle (vertex 3 is unlabelled):

```
1 100 2 101
1 100 3 -1
2 101 3 -1
```

Then, construct the `SmallGraph`:

```
SmallGraph p("pattern_graph.txt");
```

#### 2.1.2 Dynamically, using the builder methods

Construct an empty graph and add edges/anti-edges one by one:

```
SmallGraph p
  .add_edge(1, 2)
  .add_edge(1, 3)
  .add_anti_edge(2, 3)
  .set_label(1, 100)
  .set_label(2, 101)
  .set_label(3, -1);
```

### 2.2  Constructing patterns using the `PatternGenerator`

The `PatternGenerator` class is useful for

Quickly generating common patterns:
```
SmallGraph triangle = PatternGenerator::clique(3);
SmallGraph wedge = PatternGenerator::star(3);
```

Quickly generating many patterns:
```
int size = 4;
std::vector<SmallGraph> vertex_induced = PatternGenerator::all(size,
        PatternGenerator::VERTEX_BASED,        // 4 vertices
        PatternGenerator::INCLUDE_ANTI_EDGES); // anti-edges are inserted between all unadjacent vertices

std::vector<SmallGraph> vertex_induced_edge_based = PatternGenerator::all(size,
        PatternGenerator::EDGE_BASED,          // 4 edges
        PatternGenerator::INCLUDE_ANTI_EDGES); // anti-edges are inserted between all unadjacent vertices

std::vector<SmallGraph> edge_induced = PatternGenerator::all(size,
        PatternGenerator::EDGE_BASED,          // 4 edges
        PatternGenerator::EXCLUDE_ANTI_EDGES); // no extra anti-edges are inserted

std::vector<SmallGraph> edge_induced_vertex_based = PatternGenerator::all(size,
        PatternGenerator::VERTEX_BASED,        // 4 vertices
        PatternGenerator::EXCLUDE_ANTI_EDGES); // no extra anti-edges are inserted
```

Extending existing patterns:
```
std::vector<SmallGraph> vertex_extensions = PatternGenerator::extend({triangle}, PatternGenerator::VERTEX_BASED);
std::vector<SmallGraph> edge_extensions = PatternGenerator::extend({triangle}, PatternGenerator::EDGE_BASED);
```

### 2.3 Importing data graphs

Data graphs can be constructed from a preprocessed edge-list (see [Data graphs](#3-data-graphs)) or a `SmallGraph`.

```
DataGraph dg("preprocessed_data/");
```

```
SmallGraph g("small_data_graph.txt");
DataGraph dg(g);
```

### 2.4 Matching patterns

Pattern matching is done through the `match` and `count` functions.

Counting instances of patterns is very simple. Consider the following minimal motifs program:

```
#include "Peregrine.hh"
using namespace Peregrine;
void motifs(int size)
{
  int nthreads = 16;
  DataGraph g("data/citeseer/");

  auto patterns = PatternGenerator::all(size, PatternGenerator::VERTEX_BASED, PatternGenerator::INCLUDE_ANTI_EDGES);

  auto results = count(g, patterns, nthreads);

  for (const auto &[pattern, count] : results)
  {
    std::cout << pattern << ": " << count << std::endl;
  }
}
```

The arguments to `count` are straightforward: the data graph you wish to use,
the patterns whose instances you wish to count, and the number of worker
threads to use.

For arbitrary aggregations, use the `match` template function.
For example, you could replace `count` above with this snippet (note that using
`count` will be faster):

```
  const auto callback = [](auto &&handle, auto &&match) { handle.map(match.pattern, 1); };
  auto results = match<Pattern, uint64_t, AT_THE_END, UNSTOPPABLE>(g, patterns, nthreads, callback);
```

First, let's explain the template arguments to `match`:
  * `Peregrine::Pattern` is the aggregation key type
  * `uint64_t` is the aggregation value type
  * `AT_THE_END` describes when aggregation happens: either `AT_THE_END` of execution or `ON_THE_FLY`
  * `UNSTOPPABLE` indicates that you will not be using early termination

The regular parameters are the same as the `count` example except for
`callback`. This is a function that will be called on each match for a pattern.
It takes two arguments: a handle to Peregrine's aggregator, and the match itself.

Here, `callback` is mapping the pattern of the match to 1. Note that this lines
up with the types passed to `match`: the aggregation key is `Pattern` and the
value is an integer. Peregrine automatically takes care of simple types like
`uint64_t`, aggregating them with the built-in sum operator.

Using more complicated aggregation values requires them to implement a few
methods and implement a view function. The sample FSM application is a great
example: check out both `apps/fsm.cc` and `apps/Domain.hh` to see a view
function and the aggregation value methods.

Users can call three methods on the `handle`:
- `map(key, value)`: maps `key` to `value`
- `read_value(key)`: returns a view of the value mapped by `key`
- `stop()`: halts exploration early

## 3. Data graphs

Peregrine's data processor ingests graph edge-lists and stores them in binary
adjacency-list format. For labeled data graphs, a separate file of vertex-label
pairs is used. Both vertex ID's and labels should be unsigned 32-bit integers.


Edge-list file format:
```
<vertex-id> <vertex-id>
```

Label file format:
```
<vertex-id> <label>
```

Given such files `edges.txt` and `labels.txt`, the processor is used as follows:

```
$ cd Peregrine
$ mkdir data/my-graph
$ make convert_data
$ bin/convert_data edges.txt labels.txt data/my-graph/
```

To convert an unlabeled dataset, `labels.txt` can be omitted.

## 4. Reproducing our EuroSys 2020 paper results

See [the guide](experiments-guide.md).

## 5. Contributing

A sincere thank you for deciding to spend your valuable time on this project!
PR's are welcome, and will be carefully considered so long as they do not
*degrade performance* or *compromise correctness*.

If you want to understand the system, the most important bits are in:
  1. `Peregrine.hh`, which contains the main entry-points to the system
  2. `PatternMatching.hh`, which contains the pattern matching engine
  3. `Graph.hh`, where pattern analysis happens

You can also contribute in several ways without having to dig into Peregrine
internals:
- Submit your applications to be included as samples
- Point out or improve missing or unsatisfactory documentation
- Complain about confusing errors
- Suggest features you wish Peregrine had
- Report correctness or performance bugs

## 6. Acknowledgements

We are grateful for the other open-source projects this codebase relies on:
  * [bliss](http://www.tcs.hut.fi/Software/bliss/)
  * [Roaring](http://roaringbitmap.org/)

## 7. Resources

If you use this software project, please cite:

```
@inproceedings{10.1145/3342195.3387548,
  author = {Jamshidi, Kasra and Mahadasa, Rakesh and Vora, Keval},
  title = {Peregrine: A Pattern-Aware Graph Mining System},
  year = {2020},
  isbn = {9781450368827},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3342195.3387548},
  doi = {10.1145/3342195.3387548},
  booktitle = {Proceedings of the Fifteenth European Conference on Computer Systems},
  articleno = {Article 13},
  numpages = {16},
  location = {Heraklion, Greece},
  series = {EuroSys ’20}
}
```
