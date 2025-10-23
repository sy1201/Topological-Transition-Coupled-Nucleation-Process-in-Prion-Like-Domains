import numpy as np
from tqdm import tqdm
import multiprocessing

def compute_gr_vectorized(positions, box_size, r_max, dr):
    N = positions.shape[0]
    rho = N / (box_size ** 3)
    nbins = int(r_max / dr)
    r_edges = np.linspace(0, r_max, nbins + 1)
    r_centers = 0.5 * (r_edges[1:] + r_edges[:-1])

    dist_matrix = np.zeros((N, N))
    for i in range(N):
        delta = positions - positions[i]
        delta -= box_size * np.round(delta / box_size)
        dist_matrix[i] = np.linalg.norm(delta, axis=1)

    upper_triangle = dist_matrix[np.triu_indices(N, k=1)]
    hist, _ = np.histogram(upper_triangle, bins=r_edges)
    norm = (4/3) * np.pi * ((r_edges[1:]**3) - (r_edges[:-1]**3)) * rho * N
    gr = hist / norm
    return r_centers, gr

def compute_sq(r, gr, rho, q_vals):
    sq = []
    for q in q_vals:
        qr = q * r
        sin_qr_by_qr = np.sinc(qr / np.pi)
        integrand = (gr - 1) * sin_qr_by_qr * r**2
        s = 1 + 4 * np.pi * rho * np.trapz(integrand, r)
        sq.append(s)
    return np.array(sq)

def process_frame(frame_data, box_size, r_max, dr, q_vals):
    positions = np.array(frame_data) % box_size
    r, gr = compute_gr_vectorized(positions, box_size, r_max, dr)
    rho = positions.shape[0] / (box_size ** 3)
    sq = compute_sq(r, gr, rho, q_vals)
    return r, gr, q_vals, sq

def read_dump_file(filename, start_timestep, end_timestep):
    current_timestep = None
    current_atoms = []
    num_atoms = 0

    with open(filename, 'r') as f:
        line = f.readline()
        while line:
            line = line.strip()
            if "ITEM: TIMESTEP" in line:
                if current_timestep is not None and start_timestep <= current_timestep <= end_timestep:
                    if len(current_atoms) == num_atoms:
                        atoms_array = np.array(current_atoms)
                        atoms_array = atoms_array[:700 * 50]
                        chains = atoms_array.reshape((700, 50, 3))
                        centroids = np.mean(chains, axis=1)
                        yield current_timestep, centroids
                line = f.readline().strip()
                if not line: break
                current_timestep = int(line)
                current_atoms = []
                if current_timestep > end_timestep:
                    break
            elif "ITEM: NUMBER OF ATOMS" in line:
                line = f.readline().strip()
                num_atoms = int(line)
            elif "ITEM: BOX BOUNDS" in line:
                for _ in range(3):
                    f.readline()
            elif "ITEM: ATOMS" in line:
                parts = line.split()
                x_idx = parts.index('x') - 2
                y_idx = parts.index('y') - 2
                z_idx = parts.index('z') - 2
                max_idx = max(x_idx, y_idx, z_idx)
                for _ in range(num_atoms):
                    atom_line = f.readline().strip()
                    atom_data = atom_line.split()
                    x = float(atom_data[x_idx])
                    y = float(atom_data[y_idx])
                    z = float(atom_data[z_idx])
                    current_atoms.append([x, y, z])
            line = f.readline()

    if current_timestep is not None and start_timestep <= current_timestep <= end_timestep:
        if len(current_atoms) == num_atoms:
            atoms_array = np.array(current_atoms)
            atoms_array = atoms_array[:700 * 50]
            chains = atoms_array.reshape((700, 50, 3))
            centroids = np.mean(chains, axis=1)
            yield current_timestep, centroids

def main():
    box_size = 450.0
    r_max = 125.0
    dr = 5.0
    q_vals = np.linspace(0.05, 0.3, 100)
    start_timestep = 600000000
    end_timestep = 650000000

    print("正在读取dump文件...")
    num_cores = multiprocessing.cpu_count()
    parallel_jobs = max(1, num_cores // 2)
    print(f"检测到{num_cores}个CPU核心，将使用{parallel_jobs}个进程进行并行计算")

    with open("gr_700.txt", "w") as output_gr, open("sq_700.txt", "w") as output_sq:
        frame_generator = read_dump_file("dump.core.lammpstrj", start_timestep, end_timestep)
        total_frames = 0
        with open("dump.core.lammpstrj", 'r') as f:
            for line in f:
                if "ITEM: TIMESTEP" in line:
                    ts_line = f.readline().strip()
                    if ts_line.isdigit():
                        ts = int(ts_line)
                        if start_timestep <= ts <= end_timestep:
                            total_frames += 1

        if total_frames == 0:
            print("错误: 没有找到符合条件的时间步数据")
            return

        print(f"找到{total_frames}个符合条件的时间步")
        processed_frames = 0
        for ts, frame_data in tqdm(frame_generator, total=total_frames, desc="处理时间步"):
            try:
                r, gr, q, sq = process_frame(frame_data, box_size, r_max, dr, q_vals)
                output_gr.write(f"{ts}\n")
                for ri, gi in zip(r, gr):
                    output_gr.write(f"{ri:.2f} {gi:.6f}\n")
                output_sq.write(f"{ts}\n")
                for qi, si in zip(q, sq):
                    output_sq.write(f"{qi:.4f} {si:.6f}\n")
                processed_frames += 1
                del frame_data
            except Exception as e:
                print(f"处理时间步 {ts} 时出错: {e}")
                continue
    print(f"计算完成! 成功处理{processed_frames}个时间步")

if __name__ == "__main__":
    main()
