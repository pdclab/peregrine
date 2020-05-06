# How to reproduce our EuroSys'20 paper results

All these experiments were conducted on an AWS EC2 `c5.4xlarge` instance
running Ubuntu 18.04 (AMI ID: `ami-0d1cd67c26f5fca19`). To begin with, download
the
[datasets](https://drive.google.com/open?id=1GHtsboPPSDr1-nPr4kS2CfrZ5sHmEnrh)
and clone this repository.

Unfortunately, we can't provide the Friendster dataset at this time due to its
size, but the raw graph can be downloaded from
[SNAP](https://snap.stanford.edu/data/com-Friendster.html) and converted to
binary format as follows:

```
$ mkdir friendster
$ cd Peregrine
$ make convert_to_binary
$ bin/convert_to_binary.sh friendster-edges.txt 65608367 ../friendster
```

Unpack the datasets:

```
$ tar xvzf data.tgz
mico/
mico/data.bin
mico/labels.bin
orkut/
orkut/data.bin
orkut/labels.bin
patents/
patents/data.bin
patents/labels.bin
patents-labeled/
patents-labeled/data.bin
patents-labeled/labels.bin
```

The instructions below assume the following directory structure:

```
$ ls -1 .
data.tgz
friendster/
mico/
orkut/
patents/
patents-labeled/
Peregrine/
```

Next, we will build Peregrine. This requires at least GCC 9.2.1. To install
GCC 9 on Ubuntu 18.04, run:

```
$ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
$ sudo apt update
$ sudo apt install gcc-9 g++-9
$ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-9
```

Now we are ready to compile Peregrine:

```
$ cd Peregrine
$ git checkout eurosys20-experiments
$ source tbb2019/bin/tbbvars.sh intel64
$ make -j count fsm existence-query
```

Finally, we can move on to the experiments.

### Main benchmarks

We begin with the experiments that use the unchanged Peregrine code.

Note that the paper reports the mean runtimes of three executions for each
experiment; here we only show each command once.

The total execution time can be found at the end of the log file.

#### Motif counting

```
$ # 3-motifs
$ bin/count ../mico       3-motifs 16 > mi-3-motifs.log 2>&1
$ bin/count ../patents    3-motifs 16 > pa-3-motifs.log 2>&1
$ bin/count ../orkut      3-motifs 16 > ok-3-motifs.log 2>&1
$ bin/count ../friendster 3-motifs 16 > fr-3-motifs.log 2>&1
$ # 4-motifs
$ bin/count ../mico       4-motifs 16 > mi-4-motifs.log 2>&1
$ bin/count ../patents    4-motifs 16 > pa-4-motifs.log 2>&1
$ bin/count ../orkut      4-motifs 16 > ok-4-motifs.log 2>&1
```

#### FSM

```
$ # Mico
$ bin/fsm ../mico 3 2000 16 > mi-fsm-2k.log 2>&1
$ bin/fsm ../mico 3 3000 16 > mi-fsm-3k.log 2>&1
$ bin/fsm ../mico 3 4000 16 > mi-fsm-4k.log 2>&1
$ # Patents
$ bin/fsm ../patents-labeled 3 20000 16 > pa-fsm-20k.log 2>&1
$ bin/fsm ../patents-labeled 3 21000 16 > pa-fsm-21k.log 2>&1
$ bin/fsm ../patents-labeled 3 22000 16 > pa-fsm-22k.log 2>&1
$ bin/fsm ../patents-labeled 3 23000 16 > pa-fsm-23k.log 2>&1
```

#### Pattern matching

```
$ # p1
$ bin/count ../mico       query/p1.graph 16 > mi-p1.log 2>&1
$ bin/count ../patents    query/p1.graph 16 > pa-p1.log 2>&1
$ bin/count ../orkut      query/p1.graph 16 > ok-p1.log 2>&1
$ bin/count ../friendster query/p1.graph 16 > fr-p1.log 2>&1
$
$ # p2 (different labelings based on data graph)
$ bin/count ../mico       query/p2-mi.graph 16 > mi-p2.log 2>&1
$ bin/count ../patents    query/p2-pa.graph 16 > pa-p2.log 2>&1
$ bin/count ../orkut      query/p2-ok.graph 16 > ok-p2.log 2>&1
$ bin/count ../friendster query/p2-fr.graph 16 > fr-p2.log 2>&1
$
$ # p3
$ bin/count ../mico    query/p3.graph 16 > mi-p3.log 2>&1
$ bin/count ../patents query/p3.graph 16 > pa-p3.log 2>&1
$ bin/count ../orkut   query/p3.graph 16 > ok-p3.log 2>&1
$
$ # p4
$ bin/count ../mico       query/p4.graph 16 > mi-p4.log 2>&1
$ bin/count ../patents    query/p4.graph 16 > pa-p4.log 2>&1
$ bin/count ../orkut      query/p4.graph 16 > ok-p4.log 2>&1
$ bin/count ../friendster query/p4.graph 16 > fr-p4.log 2>&1
$
$ # p5
$ bin/count ../mico       query/p5.graph 16 > mi-p5.log 2>&1
$ bin/count ../patents    query/p5.graph 16 > pa-p5.log 2>&1
$ bin/count ../orkut      query/p5.graph 16 > ok-p5.log 2>&1
$ bin/count ../friendster query/p5.graph 16 > fr-p5.log 2>&1
$
$ # p6
$ bin/count ../mico       query/p6.graph 16 > mi-p6.log 2>&1
$ bin/count ../patents    query/p6.graph 16 > pa-p6.log 2>&1
$
$ # p7
$ bin/count ../mico       query/p7.graph 16 > mi-p7.log 2>&1
$ bin/count ../patents    query/p7.graph 16 > pa-p7.log 2>&1
$ bin/count ../orkut      query/p7.graph 16 > ok-p7.log 2>&1
$ bin/count ../friendster query/p7.graph 16 > fr-p7.log 2>&1
$
$ # p8
$ bin/count ../mico       query/p8.graph 16 > mi-p8.log 2>&1
$ bin/count ../patents    query/p8.graph 16 > pa-p8.log 2>&1
$ bin/count ../orkut      query/p8.graph 16 > ok-p8.log 2>&1
$ bin/count ../friendster query/p8.graph 16 > fr-p8.log 2>&1
```

#### Existence query

```
$ bin/existence-query ../mico       14-clique 16 > mi-existence.log 2>&1
$ bin/existence-query ../patents    14-clique 16 > pa-existence.log 2>&1
$ bin/existence-query ../orkut      14-clique 16 > ok-existence.log 2>&1
$ bin/existence-query ../friendster 14-clique 16 > fr-existence.log 2>&1
```

### Disabled symmetry-breaking

This experiment requires a version of the code with symmetry-breaking disabled.
Once again, the paper shows the mean of 3 executions, whereas we only present
the commands for one execution.

The log files contain the total execution time at the end.

```
$ cd Peregrine
$ make clean
$ git checkout eurosys20-prg-u
$ source tbb2019/bin/tbbvars.sh intel64
$ make -j count fsm
$
$ # 4-motifs
$ bin/count ../mico    4-motifs 16 > mi-prgu-motifs.log 2>&1
$ bin/count ../patents 4-motifs 16 > pa-prgu-motifs.log 2>&1
$ bin/count ../orkut   4-motifs 16 > ok-prgu-motifs.log 2>&1
$
$ # FSM
$ bin/fsm ../mico 3 2000 16 > mi-prgu-fsm-2k.log 2>&1
$ bin/fsm ../mico 3 3000 16 > mi-prgu-fsm-3k.log 2>&1
$ bin/fsm ../mico 3 4000 16 > mi-prgu-fsm-4k.log 2>&1
$ bin/fsm ../patents-labeled 3 19000 16 > pa-prgu-fsm-19k.log 2>&1
$ bin/fsm ../patents-labeled 3 20000 16 > pa-prgu-fsm-20k.log 2>&1
$ bin/fsm ../patents-labeled 3 21000 16 > pa-prgu-fsm-21k.log 2>&1
```

