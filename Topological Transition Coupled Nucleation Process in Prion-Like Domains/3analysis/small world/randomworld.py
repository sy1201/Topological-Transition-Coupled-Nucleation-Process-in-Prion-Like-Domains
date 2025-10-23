import numpy as np
import networkx as nx
from scipy.spatial import KDTree
import pandas as pd
import matplotlib.pyplot as plt

class LAMMPSClusterAnalyzer:
    def __init__(self, box_size=450.0, cutoff=2.0):
        self.box_size = box_size
        self.cutoff = cutoff
        self.atom_columns = {}
        self.chains = {}

    def _periodic_distance(self, x1, x2):
        dx = abs(x1 - x2)
        return min(dx, self.box_size - dx)

    def parse_timestep(self, f, ts):
        """读取一个时间步的数据"""
        self.chains.clear()
        line = f.readline()
        if not line.startswith("ITEM: NUMBER OF ATOMS"):
            raise ValueError("文件格式错误: 未找到原子数量")
        num_atoms = int(f.readline().strip())
        line = f.readline()
        if not line.startswith("ITEM: BOX BOUNDS"):
            raise ValueError("文件格式错误: 未找到边界")
        for _ in range(3):
            f.readline()
        header = f.readline().split()
        self.atom_columns = {col: idx for idx, col in enumerate(header[2:])}
        for _ in range(num_atoms):
            parts = f.readline().split()
            atom_id = int(parts[self.atom_columns['id']])
            atom_type = int(parts[self.atom_columns['type']])
            x = float(parts[self.atom_columns['x']])
            y = float(parts[self.atom_columns['y']])
            z = float(parts[self.atom_columns['z']])
            chain_id = (atom_id - 1) // 50 + 1
            if chain_id not in self.chains:
                self.chains[chain_id] = {'type1': [], 'type2': []}
            if atom_type == 1:
                self.chains[chain_id]['type1'].append((x, y, z))
            elif atom_type == 2:
                self.chains[chain_id]['type2'].append((x, y, z))

    def find_clusters(self):
        chain_ids = list(self.chains.keys())
        num_chains = len(chain_ids)
        chain_index = {cid: i for i, cid in enumerate(chain_ids)}
        parent = list(range(num_chains))

        def find(u):
            while parent[u] != u:
                parent[u] = parent[parent[u]]
                u = parent[u]
            return u

        def union(u, v):
            u_root, v_root = find(u), find(v)
            if u_root != v_root:
                parent[v_root] = u_root

        all_type1, all_type2 = [], []
        type1_chain, type2_chain = [], []
        for cid, data in self.chains.items():
            for coord in data['type1']:
                all_type1.append(coord)
                type1_chain.append(cid)
            for coord in data['type2']:
                all_type2.append(coord)
                type2_chain.append(cid)

        if len(all_type1) == 0 or len(all_type2) == 0:
            return [], []

        all_type1 = np.array(all_type1)
        all_type2 = np.array(all_type2)
        type2_tree = KDTree(all_type2)
        connections = set()

        for i, coord in enumerate(all_type1):
            indices = type2_tree.query_ball_point(coord, self.cutoff)
            min_dist, closest_idx = float('inf'), -1
            for idx in indices:
                t_coord = all_type2[idx]
                dx = self._periodic_distance(coord[0], t_coord[0])
                dy = self._periodic_distance(coord[1], t_coord[1])
                dz = self._periodic_distance(coord[2], t_coord[2])
                dist = np.sqrt(dx**2 + dy**2 + dz**2)
                if dist < min_dist:
                    min_dist, closest_idx = dist, idx
            if closest_idx != -1 and min_dist < self.cutoff:
                a, b = type1_chain[i], type2_chain[closest_idx]
                if a != b:
                    connections.add(tuple(sorted((a, b))))

        for a, b in connections:
            union(chain_index[a], chain_index[b])

        clusters = {}
        for i, cid in enumerate(chain_ids):
            root = find(i)
            if root not in clusters:
                clusters[root] = []
            clusters[root].append(cid)

        if not clusters:
            return [], []
        max_cluster = max(clusters.values(), key=len)
        return max_cluster, connections

    def build_graph(self, cluster, connections):
        cluster_set = set(cluster)
        cluster_connections = [c for c in connections if c[0] in cluster_set and c[1] in cluster_set]
        G = nx.Graph()
        G.add_nodes_from(cluster)
        G.add_edges_from(cluster_connections)
        return G

    def _component_path_length(self, G):
        if len(G) == 0:
            return 0.0
        if nx.is_connected(G):
            return nx.average_shortest_path_length(G)
        else:
            largest_cc = max(nx.connected_components(G), key=len)
            return nx.average_shortest_path_length(G.subgraph(largest_cc))

    def assess_small_world(self, G, n_rand_iter=20):
        res = {}
        n = G.number_of_nodes()
        if n == 0:
            return None
        C = nx.average_clustering(G)
        L = self._component_path_length(G)
        degrees = [d for _, d in G.degree()]
        k_mean = float(np.mean(degrees))
        m = G.number_of_edges()

        # 保度随机化
        C_rand_list, L_rand_list = [], []
        for _ in range(n_rand_iter):
            Gr = G.copy()
            nx.double_edge_swap(Gr, nswap=10*m, max_tries=100*m)
            C_rand_list.append(nx.average_clustering(Gr))
            L_rand_list.append(self._component_path_length(Gr))
        C_rand = float(np.mean(C_rand_list))
        L_rand = float(np.mean(L_rand_list))

        # 格子图
        k_even = int(round(k_mean))
        if k_even % 2 == 1:
            k_even += 1
        if k_even < 2:
            k_even = 2
        Glat = nx.watts_strogatz_graph(n, k_even, 0.0)
        C_lat = nx.average_clustering(Glat)
        L_lat = self._component_path_length(Glat)

        eps = 1e-12
        sigma = (C / max(C_rand, eps)) / (L / max(L_rand, eps))
        omega = (max(L_rand, eps) / L) - (C / max(C_lat, eps))

        res.update({
            'N': n,
            'm': m,
            'k_mean': k_mean,
            'C': C,
            'L': L,
            'C_rand': C_rand,
            'L_rand': L_rand,
            'C_lat': C_lat,
            'L_lat': L_lat,
            'C_ratio': C / max(C_rand, eps),
            'L_ratio': L / max(L_rand, eps),
            'sigma': sigma,
            'omega': omega
        })
        return res

    def run_analysis_interval(self, filename, t_start, t_end, output="testsmall.txt"):
        with open(filename, 'r') as f, open(output, 'w', encoding="utf-8") as out:
            # 英文列名
            out.write("timestep max_cluster C L C_ratio L_ratio sigma omega\n")
            while True:
                line = f.readline()
                if not line:
                    break
                if line.startswith("ITEM: TIMESTEP"):
                    ts = int(f.readline().strip())
                    if ts < t_start:
                        continue
                    if ts > t_end:
                        break
                    try:
                        self.parse_timestep(f, ts)
                        max_cluster, connections = self.find_clusters()
                        if not max_cluster:
                            continue
                        G = self.build_graph(max_cluster, connections)
                        stats = self.assess_small_world(G, n_rand_iter=30)
                        if stats is None:
                            continue
                        out.write(
                            f"{ts} {len(max_cluster)} "
                            f"{stats['C']:.6f} {stats['L']:.6f} "
                            f"{stats['C_ratio']:.4f} {stats['L_ratio']:.4f} "
                            f"{stats['sigma']:.4f} {stats['omega']:.4f}\n"
                        )
                        print(f"完成 {ts}: N={len(max_cluster)}, C={stats['C']:.6f}, L={stats['L']:.6f}, "
                              f"sigma={stats['sigma']:.4f}, omega={stats['omega']:.4f}")
                    except Exception as e:
                        print(f"跳过 {ts}: {e}")

        # 分析完成后画图
        self.plot_results(output)

    def plot_results(self, output):
        # 用 sep='\s+' 读取文件，列名去空格
        data = pd.read_csv(output, sep=r'\s+', encoding='latin1', engine='python')
        data.columns = [c.strip() for c in data.columns]

        # C/C_rand 和 L/L_rand
        plt.figure(figsize=(10,5))
        plt.plot(data["timestep"], data["C_ratio"], label="C/C_rand", marker="o")
        plt.plot(data["timestep"], data["L_ratio"], label="L/L_rand", marker="s")
        plt.axhline(1.0, color="gray", linestyle="--")
        plt.xlabel("Timestep")
        plt.ylabel("Ratio")
        plt.title("C/C_rand and L/L_rand over time")
        plt.legend()
        plt.tight_layout()
        plt.show()

        # sigma 和 最大Cluster链数 (双轴)
        fig, ax1 = plt.subplots(figsize=(10,5))
        ax1.plot(data["timestep"], data["sigma"], label="sigma", color="blue", marker="o")
        ax1.axhline(1.0, color="gray", linestyle="--")
        ax1.set_xlabel("Timestep")
        ax1.set_ylabel("sigma", color="blue")
        ax1.tick_params(axis="y", labelcolor="blue")

        ax2 = ax1.twinx()
        ax2.plot(data["timestep"], data["max_cluster"], label="Max Cluster Size", color="red", marker="s")
        ax2.set_ylabel("Max Cluster Size", color="red")
        ax2.tick_params(axis="y", labelcolor="red")

        plt.title("Sigma and Max Cluster Size over time")
        fig.tight_layout()
        plt.show()

        # omega
        plt.figure(figsize=(10,5))
        plt.plot(data["timestep"], data["omega"], label="omega", color="red", marker="o")
        plt.axhline(0.0, color="gray", linestyle="--")
        plt.xlabel("Timestep")
        plt.ylabel("omega")
        plt.title("Omega over time")
        plt.legend()
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    analyzer = LAMMPSClusterAnalyzer(box_size=450.0, cutoff=2.0)
    analyzer.run_analysis_interval(
        filename="dump.core.lammpstrj",
        t_start=600000000,
        t_end=650000000,
        output="testsmall-all.txt"
    )
