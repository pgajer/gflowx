#ifndef GRAPH_SHORTEST_PATH_H_
#define GRAPH_SHORTEST_PATH_H_

std::pair<std::vector<double>, std::vector<int>> find_all_shortest_paths_from_vertex(
    const std::vector<std::vector<int>>& adj_list,
    const std::vector<std::vector<double>>& weight_list,
    int start
	);



#endif // GRAPH_SHORTEST_PATH_H_
