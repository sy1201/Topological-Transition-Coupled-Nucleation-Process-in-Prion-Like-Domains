import numpy as np

# ----- 参数定义 -----
nlic = 50                    # 每个链亲水粒子数
ncell = 450.0
pi = np.pi
dd = 4.5                    # 相邻粒子间距
box = ncell

# ----- 计算总粒子数等 -----
np_ = nlic
mnlipid = 700
totnp = np_ * mnlipid
nr0, nc0 = 600, 600

# ----- 数据结构 -----
xnn = np.zeros((mnlipid, nlic))
ynn = np.zeros((mnlipid, nlic))
znn = np.zeros((mnlipid, nlic))

# ----- 随机链条生成 -----
for i in range(mnlipid):
    while True:
        xx, yy, zz = np.random.rand(3)
        xnn[i, 0] = xx * box
        ynn[i, 0] = yy * box
        znn[i, 0] = zz * box

        valid = True
        for j in range(1, nlic):
            v1 = np.random.rand() * pi
            v2 = np.random.rand() * 2 * pi

            xnn[i, j] = xnn[i, j - 1] + dd * np.sin(v1) * np.cos(v2)
            ynn[i, j] = ynn[i, j - 1] + dd * np.sin(v1) * np.sin(v2)
            znn[i, j] = znn[i, j - 1] + dd * np.cos(v1)

            if not (0 <= xnn[i, j] <= box and 0 <= ynn[i, j] <= box and 0 <= znn[i, j] <= box):
                valid = False
                break

        if valid:
            overlap = False
            for j in range(nlic - 1):
                dx = xnn[i, j + 1] - xnn[i, j]
                dy = ynn[i, j + 1] - ynn[i, j]
                dz = znn[i, j + 1] - znn[i, j]
                rr = dx ** 2 + dy ** 2 + dz ** 2
                if rr < 0.1:
                    overlap = True
                    break
            if not overlap:
                break

# ----- 输出 LAMMPS 数据文件 -----
bonds = (nlic - 1) * mnlipid

with open("fus.data", "w") as f:
    f.write("data file\n\n")
    f.write(f"{mnlipid * nlic:6d} atoms\n")
    f.write(f"{bonds:6d} bonds\n")
    f.write(f"{0:6d} angles\n")
    f.write(f"{0:6d} dihedrals\n")
    f.write(f"{0:6d} impropers\n\n")

    f.write(f"{3:6d} atom types\n")
    f.write(f"{1:6d} bond types\n")
    f.write(f"{0:6d} angle types\n")
    f.write(f"{0:6d} dihedral types\n")
    f.write(f"{0:6d} improper types\n\n")

    f.write(f"{0.0:7.2f} {ncell:7.2f} xlo xhi\n")
    f.write(f"{0.0:7.2f} {ncell:7.2f} ylo yhi\n")
    f.write(f"{0.0:7.2f} {ncell:7.2f} zlo zhi\n\n")

    f.write("Masses\n\n")
    for i in range(1, 4):
        f.write(f"{i} 3534.291740\n")

    f.write("\nAtoms\n\n")
    for i in range(mnlipid):
        for j in range(7, 30, 8):
            k = i * np_ + j + 1  # 粒子ID从1开始
            f.write(f"{k:6d} {i+1:6d} {1:6d} {xnn[i,j]:12.4f} {ynn[i,j]:12.4f} {znn[i,j]:12.4f}\n")
        for j in range(31, nlic, 8):
            k = i * np_ + j + 1
            f.write(f"{k:6d} {i+1:6d} {2:6d} {xnn[i,j]:12.4f} {ynn[i,j]:12.4f} {znn[i,j]:12.4f}\n")
        for m in range(1, 8):
            for j in range(m - 1, nlic, 8):
                k = i * np_ + j + 1
                f.write(f"{k:6d} {i+1:6d} {3:6d} {xnn[i,j]:12.4f} {ynn[i,j]:12.4f} {znn[i,j]:12.4f}\n")

    f.write("\nBonds\n\n")
    bond_id = 1
    for m in range(mnlipid):
        for i in range(nlic - 1):
            atom1 = m * np_ + i + 1
            atom2 = atom1 + 1
            f.write(f"{bond_id:6d} {1:7d} {atom1:7d} {atom2:7d}\n")
            bond_id += 1
