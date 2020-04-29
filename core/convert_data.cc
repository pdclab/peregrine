#include <unistd.h>
#include <assert.h>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "tbb/parallel_sort.h"
#include <thread>
#include <sys/mman.h>
#include <condition_variable>
#include <mutex>
#include <sys/time.h>
#include <iostream>
#include <atomic>
#include <chrono>
#include <thread>

struct adjlist {
  adjlist() : length(0), ptr(nullptr) {}
  adjlist(uint32_t l, uint32_t *p) : length(l), ptr(p) {}
  std::atomic<uint32_t> length;
  uint32_t *ptr;
};

bool is_directory(const std::string &path)
{
   struct stat statbuf;
   if (stat(path.c_str(), &statbuf) != 0)
       return 0;
   return S_ISDIR(statbuf.st_mode);
}


uint64_t max_thread_count = std::thread::hardware_concurrency();
std::string delimiters = " ,\t";
uint64_t MAX_V_ID=130'000'000;
std::vector<std::vector<uint32_t>> per_thread_degree_map(max_thread_count);
std::vector<uint32_t *> per_thread_edge_list(max_thread_count);
std::vector<uint32_t> per_thread_edge_cursor(max_thread_count,0);
std::vector<std::pair<uint32_t,uint32_t>> degree_map;
std::vector<uint32_t> all_degree;
std::vector<uint32_t> vertex_adj_map;
std::vector<uint32_t> vertex_degree_map;
std::vector<uint32_t> vertex_start_offset;
adjlist *data_graph;
uint32_t vertex_id=1;
size_t bytes_written;

template<typename T>
void empty_swap(std::vector<T>& vec) {
   std::vector<T>().swap(vec);
}

void calculate_degree_map(uint32_t thread_id, uint64_t file_block_size, char* graph_data, uint32_t num_vertices){
  per_thread_degree_map[thread_id].resize(num_vertices,0);
  uint64_t cursor_start = thread_id*file_block_size;
  uint64_t cursor_end = (thread_id+1)*file_block_size;
  if(thread_id!=0){
    while(graph_data[cursor_start] != '\n') cursor_start++;
    cursor_start++;
  }
  if(thread_id!=max_thread_count-1){
    while(graph_data[cursor_end] != '\n') cursor_end++;
  }
  while (cursor_start < cursor_end){
    std::string a = "";
    std::string b = "";
    uint32_t s1, s2;
    while (graph_data[cursor_start] == '#'){
      while (graph_data[cursor_start] != '\n')
      cursor_start++;
      cursor_start++;
    }
    while (delimiters.find(graph_data[cursor_start]) == std::string::npos){
      a += graph_data[cursor_start];
      cursor_start++;
    }
    while (delimiters.find(graph_data[cursor_start]) != std::string::npos)
      cursor_start++;
    while (graph_data[cursor_start] != '\n'){
      b += graph_data[cursor_start];
      cursor_start++;
    }
    cursor_start++;
    s1 = std::atoi(a.c_str());
    s2 = std::atoi(b.c_str());
    if (s1 > s2 ){
      uint32_t temp = s1;
      s1 = s2;
      s2 = temp;
    }

    if (s2 >= MAX_V_ID)
    {
      std::cerr << "ERROR: Observed a vertex in the edge list with ID " << s2 << " > maximum vertex ID " << MAX_V_ID-1 << "." << std::endl;
      std::cerr << "Exiting..." << std::endl;
      exit(-1);
    }


    if(s1!=s2){
      per_thread_degree_map[thread_id][s1]++;
      per_thread_degree_map[thread_id][s2]++;
    }
  }
}

void populate_data_graph(uint32_t thread_id, uint64_t file_block_size, char* graph_data){

  uint64_t cursor_start = thread_id*file_block_size;
  uint64_t cursor_end = (thread_id+1)*file_block_size;
  if(thread_id!=0){
    while(graph_data[cursor_start] != '\n') cursor_start++;
    cursor_start++;
  }
  if(thread_id!=max_thread_count-1){
    while(graph_data[cursor_end] != '\n') cursor_end++;
  }
  if(thread_id==0){
    while (graph_data[cursor_start] == '#'){
      while (graph_data[cursor_start] != '\n')
        cursor_start++;
      cursor_start++;
    }
  }
  while (cursor_start < cursor_end){
    std::string a = "";
    std::string b = "";
    uint32_t s1, s2;
    //skip comments
    while (delimiters.find(graph_data[cursor_start]) == std::string::npos){
      a += graph_data[cursor_start];
      cursor_start++;
    }
    while (delimiters.find(graph_data[cursor_start]) != std::string::npos)
      cursor_start++;
    while (graph_data[cursor_start] != '\n'){
      b += graph_data[cursor_start];
      cursor_start++;
    }
    cursor_start++;
    s1 = std::atoi(a.c_str());
    s2 = std::atoi(b.c_str());
    if(s1!=s2){
      uint32_t curr_pos_s2 = data_graph[vertex_adj_map[s2]].length++;
      data_graph[vertex_adj_map[s2]].ptr[curr_pos_s2] = s1;
      uint32_t curr_pos_s1 = data_graph[vertex_adj_map[s1]].length++;
      data_graph[vertex_adj_map[s1]].ptr[curr_pos_s1] = s2;
    }
    else{
        std::cout<<"Loop Found "<<s1<<std::endl;
    }
  }
}

void write_graph_to_disk(uint32_t start, uint32_t end, int adj_file,std::atomic<uint32_t> &edge_count,std::atomic<uint32_t> &data_count){
    uint64_t curr_vertex_offset = 0;
    if(start == 0) curr_vertex_offset = vertex_start_offset[0]*4;
    uint64_t max_buffer = 3*256*512*512;
    uint32_t *stage = new uint32_t[max_buffer+1];
    uint32_t stage_cursor = 0;
    uint32_t local_data_count=0;
    uint32_t local_edge_count=0;
    for(uint32_t i = start; i < end; i++){
      uint32_t curr_v = vertex_adj_map[degree_map[i].first];
      uint32_t source = vertex_degree_map[degree_map[i].first];
      uint32_t degree = degree_map[i].second;
      uint32_t new_degree = degree;
      local_edge_count+=new_degree;

      stage[stage_cursor] = source;
      stage_cursor++;

      local_data_count++;

      stage[stage_cursor] = new_degree;
      stage_cursor++;

      local_data_count++;

      for(uint32_t j = 0; j < new_degree;j++){
        local_data_count++;
        stage[stage_cursor] = data_graph[curr_v].ptr[j];
        stage_cursor++;
      }

      if(stage_cursor > 512*512*512){
        bytes_written = pwrite(adj_file,stage,sizeof(uint32_t)*(stage_cursor),curr_vertex_offset);
        curr_vertex_offset += sizeof(uint32_t)*(stage_cursor);
        stage_cursor = 0;
      }
    }

    if(stage_cursor != 0){
      bytes_written = pwrite(adj_file,stage,sizeof(uint32_t)*(stage_cursor),curr_vertex_offset);
    }
    data_count += local_data_count;
    edge_count += local_edge_count;

}

int main(int argc, char* argv[]){
  if(argc < 4){
    std::cout<<"Usage : " << argv[0] << " <input edge file> [input label file] <max vertex ID> <output directory>"<<std::endl;
    exit(0);
  }

  std::string data_graph_path(argv[1]);
  std::string outdir;
  std::string label_path;
  if (argc == 4)
  {
    // no labels
    MAX_V_ID = std::atoi(argv[2]) + 1;
    outdir = argv[3];
  }
  else
  {
    label_path = argv[2];
    MAX_V_ID = std::atoi(argv[3]) + 1;
    outdir = argv[4];
  }

  if (!is_directory(outdir))
  {
    std::cerr << "ERROR: " << outdir << " is not a directory." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    exit(-1);
  }

  // check for MAX_V_ID overflow
  if (MAX_V_ID-1 >= static_cast<uint32_t>(-1))
  {
    std::cerr << "ERROR: Provided maximum vertex ID overflows." << std::endl;
    std::cerr << "Exiting..." << std::endl;
    exit(-1);
  }

  char* graph_data;

  std::chrono::steady_clock::time_point begin,end,g_begin,g_end;
  std::thread workers[max_thread_count];
  g_begin = std::chrono::steady_clock::now();

  struct stat st;

  begin = std::chrono::steady_clock::now();

  stat(data_graph_path.c_str(), &st);
  size_t input_file_size = st.st_size;
  std::cout<<"Size of file is "<<input_file_size<<std::endl;
  int fd = open(data_graph_path.c_str(), O_RDONLY, 0);
  assert(fd != -1);

  graph_data = static_cast<char*>(mmap(NULL, input_file_size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0));
  assert(graph_data != MAP_FAILED);
  end = std::chrono::steady_clock::now();
    std::cout << "Time Taken for memory mapping the file : " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
  uint64_t file_block_size = input_file_size/max_thread_count;
  begin = std::chrono::steady_clock::now();

  for(int i = 0; i < max_thread_count; i++) {
    workers[i] = std::thread(calculate_degree_map,i,file_block_size,graph_data, MAX_V_ID);
  }

  for(int i = 0; i < max_thread_count; i++) workers[i].join();

  uint64_t total_degree=0;
  all_degree.resize(MAX_V_ID,0);
  vertex_degree_map.resize(MAX_V_ID,0);
  vertex_adj_map.resize(MAX_V_ID,0);
  std::cout<<"degree map calculated per thread"<<std::endl;
  end = std::chrono::steady_clock::now();
  std::cout << "Time Taken degree calculation per thread " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
  begin = std::chrono::steady_clock::now();

  for(uint32_t thread_id = 0; thread_id < max_thread_count;thread_id++){
    for(uint32_t j = 0; j < per_thread_degree_map[thread_id].size(); j++){
      uint32_t curr_source = j;
      uint32_t curr_degree = per_thread_degree_map[thread_id][j];
      total_degree += curr_degree;
      all_degree[j] += curr_degree;
    }
    empty_swap(per_thread_degree_map[thread_id]);
  }

  for(uint32_t i = 0; i < MAX_V_ID;i++){
    if(all_degree[i]!=0){
      vertex_degree_map[i] = vertex_id;
      vertex_adj_map[i] = vertex_id-1;
      degree_map.push_back({i,all_degree[i]});
      vertex_id++;
    }
  }
  empty_swap(all_degree);

  end = std::chrono::steady_clock::now();
  std::cout << "Time Taken global degree calculation  " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;


  std::cout<<"degree map calculated"<<std::endl;
  std::cout<<"Number of vertices : "<<vertex_id-1<<std::endl;
  std::cout<<"Total Degree : "<<total_degree<<std::endl;
  begin = std::chrono::steady_clock::now();

  // for(int i = 0; i < max_thread_count; i++){
  //     workers[i] = std::thread(populate_data_graph,i,file_block_size);
  // }
  // for(int i = 0; i < max_thread_count; i++) workers[i].join();




  data_graph = new adjlist[vertex_id-1];
  tbb::parallel_for( tbb::blocked_range<int>(0,vertex_id-1),
             [&](tbb::blocked_range<int> r){
    for(int i = r.begin(); i < r.end(); ++i){
      data_graph[i].ptr = new uint32_t[degree_map[i].second];
      data_graph[i].length = 0;
    }
  });

  for(int i=0;i<max_thread_count;i++){
    workers[i] = std::thread(populate_data_graph,i,file_block_size,graph_data);
  }
  for(int i=0;i<max_thread_count;i++) workers[i].join();


  munmap(graph_data, input_file_size);
  close(fd);

  std::cout<<"data graph  created"<<std::endl;



  end = std::chrono::steady_clock::now();
  std::cout << "Time Taken for populating data  " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
  std::cout<<"data graph populated"<<std::endl;
  tbb::parallel_sort(degree_map.begin(),degree_map.end(),
  [=]( const std::pair< uint32_t, uint32_t >&a, const std::pair< uint32_t, uint32_t >&b ){
    if (a.second == b.second) return a.first < b.first;
    return a.second > b.second;
  });
  std::cout<<"degree map sorted"<<std::endl;

  tbb::parallel_for( tbb::blocked_range<int>(0,vertex_id-1),
             [&](tbb::blocked_range<int> r){
    for(int i = r.begin(); i < r.end(); ++i){
      vertex_degree_map[degree_map[i].first] = i+1;
    }
  });

  std::cout<<"Biggest Degree of the vertex : "<<degree_map[0].first<<" "<<degree_map[0].second<<std::endl;
  uint64_t offset=0;
  std::atomic<uint32_t> data_count{0};
  FILE* adj_file[max_thread_count];

  for(uint32_t i = 0; i < max_thread_count; i++)
    adj_file[i] = fopen((outdir+"/data_"+std::to_string(i)+".bin").c_str(),"wb+");
  vertex_id--;
  std::atomic<uint32_t> edge_count{0};
  uint32_t vertex_start_idx = 3+vertex_id;
  bytes_written = pwrite(fileno(adj_file[0]),&vertex_start_idx,sizeof(uint32_t),offset); // place holder for total file size
  offset+=sizeof(uint32_t);
  data_count++;
  bytes_written = pwrite(fileno(adj_file[0]),&vertex_id,sizeof(uint32_t),offset); // write number of vertices
  offset+=sizeof(uint32_t);
  data_count++;
  bytes_written = pwrite(fileno(adj_file[0]),&vertex_start_idx,sizeof(uint32_t),offset); // place holder for edge count
  offset+=sizeof(uint32_t);
  data_count++;

  vertex_start_offset.resize(vertex_id,0);
  tbb::parallel_for( tbb::blocked_range<int>(0,vertex_id),
            [&](tbb::blocked_range<int> r){
    uint32_t start = r.begin();
    uint32_t end = r.end();
    for(uint32_t i = start; i < end; i++){
    uint32_t curr_v = vertex_adj_map[degree_map[i].first];
    uint32_t degree = degree_map[i].second;
    assert(degree == data_graph[curr_v].length);
    std::vector<uint32_t> curr_neighbours;
    curr_neighbours.resize(degree);
    for(uint32_t k=0; k < degree; k++) curr_neighbours[k] = vertex_degree_map[data_graph[curr_v].ptr[k]];
    std::sort(curr_neighbours.begin(),curr_neighbours.end());
    auto last = std::unique(curr_neighbours.begin(), curr_neighbours.end());
    curr_neighbours.erase(last,curr_neighbours.end());
    uint32_t new_degree=curr_neighbours.size();
    for(uint32_t j = 0; j < new_degree; j++){
      data_graph[curr_v].ptr[j] = curr_neighbours[j];
    }
    data_graph[curr_v].length = new_degree;
    degree_map[i].second = new_degree;
    }
  });

  uint32_t* v_offset = new uint32_t[vertex_id];
  vertex_start_offset[0] = 3+vertex_id;
  v_offset[0] = 3+vertex_id;
  data_count++;

  for(uint32_t i=1; i < vertex_id; i++){
    vertex_start_offset[i] = vertex_start_offset[i-1]+2+degree_map[i-1].second;
    v_offset[i] = vertex_start_offset[i];
    data_count++;
  }
  bytes_written = pwrite(fileno(adj_file[0]),v_offset,sizeof(uint32_t)*vertex_id,offset);
  offset += sizeof(uint32_t)*vertex_id;
  std::vector<std::pair<uint32_t,uint32_t>> v_work;
  uint32_t v_start;
  uint32_t uniform_e = total_degree/max_thread_count;
  v_start = 0;
  uint32_t curr_e_count = 0;
  for(uint32_t i = 0; i < vertex_id; i++){
    curr_e_count += degree_map[i].second;
    if(curr_e_count >= uniform_e){
      curr_e_count = 0;
      if(v_work.size() == max_thread_count-1){
          v_work.push_back({v_start,vertex_id});
          break;
      }
      v_work.push_back({v_start,i+1});
      v_start = i+1;
    }
    if(i == vertex_id - 1){
      v_work.push_back({v_start,vertex_id});
      break;
    }
  }
  begin = std::chrono::steady_clock::now();
  for(uint32_t i = 0; i < max_thread_count; i++){
      workers[i] = std::thread(write_graph_to_disk,v_work[i].first,v_work[i].second,fileno(adj_file[i]),std::ref(edge_count),std::ref(data_count));
  }
  for(uint32_t i = 0; i < max_thread_count; i++) workers[i].join();
  end = std::chrono::steady_clock::now();
  std::cout << "Time Taken for writing data to disk  " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

  uint32_t stage_value = data_count.load();
  bytes_written = pwrite(fileno(adj_file[0]),&stage_value,sizeof(uint32_t),0); // fill the gap
  stage_value = edge_count.load();
  stage_value/=2;
  bytes_written = pwrite(fileno(adj_file[0]),&stage_value,sizeof(uint32_t),8); // fill the gap
  std::cout<<"Edge Count : "<<edge_count<<" "<<"Data Count : "<<data_count<<std::endl;
  for(uint32_t i = 0; i < max_thread_count; i++) close(fileno(adj_file[i]));


  if (argc > 4) {

    std::cout<<"Preparing Labels data"<<std::endl;
    struct stat st;
    stat(label_path.c_str(), &st);
    size_t labels_file_size = st.st_size;
    std::cout<<"Size of Label file is "<<labels_file_size<<std::endl;
    int labels_fd = open(label_path.c_str(), O_RDONLY, 0);
    assert(labels_fd != -1);
    FILE* labels_file = fopen((outdir+"/labels.bin").c_str(),"wb+");
    char* labels_data = static_cast<char*>(mmap(NULL, labels_file_size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, labels_fd, 0));
    off_t labels_offset=0;
    file_block_size = labels_file_size/max_thread_count;
    for(uint32_t i = 0; i < max_thread_count;i++){
      uint64_t cursor_start = i*file_block_size;
      uint64_t cursor_end = (i+1)*file_block_size;
      if(i!=0){
        while(labels_data[cursor_start] != '\n') cursor_start++;
        cursor_start++;
      }
      if(i!=max_thread_count-1){
        while(labels_data[cursor_end] != '\n') cursor_end++;
      }
      while (cursor_start < cursor_end)
      {
        std::string a = "";
        std::string b = "";
        uint32_t s1, s2;
        //skip comments
        while (labels_data[cursor_start] == '#')
        {
          while (labels_data[cursor_start] != '\n')
          cursor_start++;
          cursor_start++;
        }
        while (delimiters.find(labels_data[cursor_start]) == std::string::npos)
        {
          a += labels_data[cursor_start];
          cursor_start++;
        }
        while (delimiters.find(labels_data[cursor_start]) != std::string::npos) cursor_start++;
        while (labels_data[cursor_start] != '\n')
        {
          b += labels_data[cursor_start];
          cursor_start++;
        }
        cursor_start++;
        s1 = std::atoi(a.c_str());
        s2 = std::atoi(b.c_str());
        uint32_t curr_vertex = vertex_degree_map[s1];
        if(curr_vertex!=0){
	bytes_written = pwrite(fileno(labels_file),&curr_vertex,sizeof(uint32_t),labels_offset);
        labels_offset+=sizeof(uint32_t);
        bytes_written = pwrite(fileno(labels_file),&s2,sizeof(uint32_t),labels_offset);
        labels_offset+=sizeof(uint32_t);
}
      }
    }
    munmap(labels_data, labels_file_size);
    close(labels_fd);
    close(fileno(labels_file));

  }
g_end = std::chrono::steady_clock::now();
  std::cout << "Total Time Taken : " << std::chrono::duration_cast<std::chrono::microseconds>(g_end - g_begin).count() << "[µs]" << std::endl;

  return 0;
}
