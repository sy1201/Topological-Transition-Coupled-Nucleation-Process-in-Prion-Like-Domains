import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from tqdm import tqdm
from scipy.spatial import KDTree
import os

# 设置中文字体
plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]
plt.rcParams["axes.unicode_minus"] = False

class LAMMPSClusterAnalyzer:
    def __init__(self, box_size=450.0, cutoff=2.0, time_step=629300000):
        self.box_size = box_size
        self.cutoff = cutoff
        self.chains = {}
        self.time_step = time_step

    def read_lammpstrj(self, filename):
        print(f"正在读取文件: {filename}")
        with open(filename, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        timestep_pos = -1
        for i, line in enumerate(lines):
            if line.startswith(f"{self.time_step}"):
                timestep_pos = i
                break
        if timestep_pos == -1:
            raise ValueError(f"文件中未找到时间步 {self.time_step}")
        
        if lines[timestep_pos + 1] != "ITEM: NUMBER OF ATOMS":
            raise ValueError("文件格式错误: 未找到原子数量")
        num_atoms = int(lines[timestep_pos + 2])
        
        atoms_start = -1
        for i in range(timestep_pos + 3, len(lines)):
            if lines[i].startswith("ITEM: ATOMS"):
                atoms_header = lines[i].split()
                self.atom_columns = {col: idx for idx, col in enumerate(atoms_header[2:])}
                atoms_start = i + 1
                break
        if atoms_start == -1:
            raise ValueError("文件格式错误: 未找到原子数据")
        
        atoms_end = atoms_start + num_atoms
        required_columns = ['id', 'type', 'x', 'y', 'z']
        for col in required_columns:
            if col not in self.atom_columns:
                raise ValueError(f"缺少必要列: {col}")
        
        type_count = {1:0, 2:0, 3:0}
        for line in tqdm(lines[atoms_start:atoms_end], desc="解析原子数据"):
            parts = line.split()
            if len(parts) < len(self.atom_columns):
                continue
            atom_id = int(parts[self.atom_columns['id']])
            atom_type = int(parts[self.atom_columns['type']])
            x = float(parts[self.atom_columns['x']])
            y = float(parts[self.atom_columns['y']])
            z = float(parts[self.atom_columns['z']])
            if atom_type in type_count:
                type_count[atom_type] += 1
            chain_id = (atom_id - 1) // 50 + 1
            if chain_id not in self.chains:
                self.chains[chain_id] = {'type1': [], 'type2': [], 'type3': []}
            if atom_type == 1:
                self.chains[chain_id]['type1'].append((x,y,z))
            elif atom_type == 2:
                self.chains[chain_id]['type2'].append((x,y,z))
            elif atom_type == 3:
                self.chains[chain_id]['type3'].append((x,y,z))
        print(f"解析完成，共 {len(self.chains)} 条链")
        print(f"类型1: {type_count[1]}, 类型2: {type_count[2]}, 类型3: {type_count[3]}")
        if type_count[1] == 0 or type_count[2] == 0:
            raise ValueError("缺少类型1或类型2的珠子，无法分析cluster")

    def _periodic_distance(self, x1, x2):
        dx = abs(x1 - x2)
        return min(dx, self.box_size - dx)

    def find_clusters(self):
        chain_ids = list(self.chains.keys())
        num_chains = len(chain_ids)
        chain_index = {cid:i for i,cid in enumerate(chain_ids)}
        parent = list(range(num_chains))
        
        def find(u):
            while parent[u] != u:
                parent[u] = parent[parent[u]]
                u = parent[u]
            return u
        def union(u,v):
            u_root = find(u)
            v_root = find(v)
            if u_root != v_root:
                parent[v_root] = u_root

        # 准备类型1和类型2粒子
        all_type1, all_type2 = [], []
        type1_chain, type2_chain = [], []
        for cid, data in self.chains.items():
            for coord in data['type1']:
                all_type1.append(coord)
                type1_chain.append(cid)
            for coord in data['type2']:
                all_type2.append(coord)
                type2_chain.append(cid)
        all_type1 = np.array(all_type1)
        all_type2 = np.array(all_type2)
        type2_tree = KDTree(all_type2)
        connections = set()
        print("正在计算链条间连接（每条链一对一最短距离）...")
        for i, coord in enumerate(tqdm(all_type1, desc="type1粒子")):
            indices = type2_tree.query_ball_point(coord, self.cutoff)
            min_dist = float('inf')
            closest_idx = -1
            for idx in indices:
                t_coord = all_type2[idx]
                dx = self._periodic_distance(coord[0], t_coord[0])
                dy = self._periodic_distance(coord[1], t_coord[1])
                dz = self._periodic_distance(coord[2], t_coord[2])
                dist = np.sqrt(dx**2 + dy**2 + dz**2)
                if dist < min_dist:
                    min_dist = dist
                    closest_idx = idx
            if closest_idx != -1 and min_dist < self.cutoff:
                a = type1_chain[i]
                b = type2_chain[closest_idx]
                if a != b:
                    connections.add(tuple(sorted((a,b))))

        # Union-Find 构建cluster
        for a,b in tqdm(connections, desc="构建Cluster"):
            union(chain_index[a], chain_index[b])

        clusters = {}
        for i, cid in enumerate(chain_ids):
            root = find(i)
            if root not in clusters:
                clusters[root] = []
            clusters[root].append(cid)
        max_cluster = max(clusters.values(), key=len)
        print(f"最大Cluster包含 {len(max_cluster)} 条链")
        return max_cluster, connections

    def build_graph(self, cluster, connections):
        cluster_set = set(cluster)
        cluster_connections = [c for c in connections if c[0] in cluster_set and c[1] in cluster_set]
        G = nx.Graph()
        G.add_nodes_from(cluster)
        G.add_edges_from(cluster_connections)
        return G

    def calculate_betweenness(self, G):
        print("计算介数中心性...")
        return nx.betweenness_centrality(G, normalized=True), nx.betweenness_centrality(G, normalized=False)

    def calculate_small_world_properties(self, G):
        print("计算小世界网络特性...")
        C = nx.average_clustering(G)
        if nx.is_connected(G):
            L = nx.average_shortest_path_length(G)
        else:
            largest_cc = max(nx.connected_components(G), key=len)
            L = nx.average_shortest_path_length(G.subgraph(largest_cc))
            print(f"图不连通，使用最大连通组件计算平均路径长度 ({len(largest_cc)} 节点)")
        return C, L

    def save_results(self, max_cluster, G, betweenness_norm, betweenness_raw):
        # 保存最大Cluster信息
        with open('maxtime-140.txt', 'w', encoding='utf-8') as f:
            f.write(f"最大Cluster链数: {len(max_cluster)}\n")
            f.write("链序号: " + ", ".join(map(str, sorted(max_cluster))) + "\n")
        
        # 保存理论数据，包括归一化和未归一化的介数中心性
        with open('theory-140.txt', 'w', encoding='utf-8') as f:
            f.write("节点 度数 未归一化介数中心性 归一化介数中心性\n")
            for n in sorted(G.nodes()):
                f.write(f"{n} {G.degree(n)} {betweenness_raw[n]:.6f} {betweenness_norm[n]:.6f}\n")
        
        # 度数为6的节点介数中心性分布（使用未归一化值，可保留原逻辑）
        dist_values = [betweenness_raw[n] for n in G.nodes() if G.degree(n) == 5]
        if dist_values:
            bin_width = 50
            max_val = max(dist_values)
            bins = np.arange(0, max_val + bin_width, bin_width)
            hist, edges = np.histogram(dist_values, bins=bins)
            prob = hist / hist.sum()
            with open('distribution-140.txt', 'w', encoding='utf-8') as f:
                f.write("介数中心性\t概率\n")
                for left, p in zip(edges[:-1], prob):
                    f.write(f"{left:.2f}\t{p:.6f}\n")
            # 绘制直方图
            plt.figure(figsize=(8,5))
            plt.bar(edges[:-1], prob, width=bin_width, align='edge', color='skyblue', edgecolor='black')
            plt.xlabel("介数中心性")
            plt.ylabel("概率")
            plt.title("度数为6节点介数中心性分布")
            plt.tight_layout()
            plt.savefig("distribution_hist.png", dpi=300)
            plt.show()
        else:
            with open('distribution-140.txt', 'w', encoding='utf-8') as f:
                f.write("没有度数为6节点\n")

    def save_small_world_results(self, C, L, num_nodes):
        with open('small world-140.txt', 'w', encoding='utf-8') as f:
            f.write("小世界网络特性分析\n")
            f.write(f"节点数: {num_nodes}\n")
            f.write(f"平均聚集系数 C: {C:.6f}\n")
            f.write(f"平均路径长度 L: {L:.6f}\n")
            f.write("参考: 高聚集系数且平均路径长度短\n")
        print("小世界网络结果已保存")

    def visualize_graph(self, G, betweenness_norm):
        print("绘制网络图...")
        num_nodes = len(G.nodes())
        pos = nx.spring_layout(G, seed=42, k=0.8, iterations=800)
        nodesize = 50

        # 使用归一化的介数中心性，并限制最大值为0.2
        bc_values = np.array(list(betweenness_norm.values()))
        bc_values = np.clip(bc_values, 0, 0.2)

        norm = Normalize(vmin=0, vmax=0.2)
        cmap = LinearSegmentedColormap.from_list('blue_to_yellow', ['blue', 'cyan', 'green', 'yellow'])
        node_colors = cmap(norm(bc_values))

        fig, ax = plt.subplots(figsize=(8,6))
        nx.draw_networkx_edges(G, pos, edge_color='lightskyblue', width=1.0, alpha=0.7, ax=ax)
        nx.draw_networkx_nodes(G, pos, node_size=nodesize, node_color=node_colors, alpha=1.0,
                            edgecolors='black', linewidths=0.5, ax=ax)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)

        # 修改 colorbar 标签字体
        cbar.set_label(
            'Normalized Betweenness Centrality',
            fontsize=22,              # 标签字号
            fontweight='bold',        # 加粗
            family='Times New Roman'  # 字体
        )

        # 修改刻度字体
        cbar.ax.tick_params(labelsize=16)  # 刻度字号
        for tick in cbar.ax.get_yticklabels():
            tick.set_fontname('Times New Roman')
            tick.set_fontweight('bold')

        plt.axis('off')
        plt.tight_layout()
        plt.savefig('cluster_network-100.png', dpi=300)
        plt.show()


    def run_analysis(self, filename):
        try:
            self.read_lammpstrj(filename)
            max_cluster, connections = self.find_clusters()
            G = self.build_graph(max_cluster, connections)
            betweenness_norm, betweenness_raw = self.calculate_betweenness(G)
            self.save_results(max_cluster, G, betweenness_norm, betweenness_raw)
            C, L = self.calculate_small_world_properties(G)
            self.save_small_world_results(C, L, len(G.nodes()))
            self.visualize_graph(G, betweenness_norm)
            print("分析完成！")
        except Exception as e:
            print(f"分析过程中出错: {str(e)}")
            raise

if __name__ == "__main__":
    analyzer = LAMMPSClusterAnalyzer(box_size=450.0, cutoff=2.0)
    input_file = "dump.core.lammpstrj"
    if os.path.exists(input_file):
        analyzer.run_analysis(input_file)
    else:
        print(f"错误: 找不到文件 {input_file}")
