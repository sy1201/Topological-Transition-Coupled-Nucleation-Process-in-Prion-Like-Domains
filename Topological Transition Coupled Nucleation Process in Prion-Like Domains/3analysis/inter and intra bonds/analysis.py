import datetime
import dumpreader as dr
import argparse
import math_function as mf
import numpy as np
import collections

parser = argparse.ArgumentParser()
parser.add_argument("lammpstrj", help="input full name of .lammpstrj file!")
parser.add_argument("-bm",'--box_mid', help="let box center [0,0,0]",action="store_true")
parser.add_argument("-cm",'--cluster_mid', help="let max cluster center of mass [0,0,0]",action="store_true")
parser.add_argument("-yzs",'--cluster_yzs', help="use yzs_system",action="store_true")
parser.add_argument("-sy",'--cluster_sy', help="use sy_system",action="store_true")
#analysis
parser.add_argument("-rg",'--mol_rg', help="output every mol rg",action="store_true")
parser.add_argument("-c",'--cluster', help="output cluster message",action="store_true")
parser.add_argument("-mc",'--max_cluster', help="output the max_cluster message",action="store_true")
parser.add_argument("-mcb",'--max_interbonds', help="output the max_cluster bonds message",action="store_true")
parser.add_argument("-cv",'--cluster_volum', help="output the cluster zise statistics results",action="store_true")


parser.add_argument("-o",'--output', help="output new lammpstrj",type=str)
parser.add_argument("-r",'--read', help="read freame from start to end every step",action="extend", nargs="+", type=int)



args = parser.parse_args()
#配置-r
frame_strat = 0
frame_end = -1
frame_step = 1
if not args.read:
    print('note read all frame!')
elif len(args.read) == 3:
    frame_strat = args.read[0]
    frame_end = args.read[1]
    frame_step = args.read[2]
elif  len(args.read) == 2:
    frame_strat = args.read[0]
    frame_end = args.read[1]
elif  len(args.read) == 1:
    frame_end = args.read[0]
else:
    exit('wrong -r commend!')

if args.cluster_yzs and args.cluster_sy:
    exit('choose only 1 mode!-yzs or -sy')
if not (args.cluster_yzs or args.cluster_sy):
    exit('choose -yzs or -sy for cluster')

cut_off = 3.5
num_chain = 700 # FUS 链条数
num_beads = 50 # 每个FUS 珠子数




def main():

    start=datetime.datetime.now()

    # mass1=[-1,131.199997,71.0800018,87.0800018,114.099998,115.099998,163.199997,101.099998,128.100006,57.0499992,97.1200027]
    # open file 初始化体系
    r = dr.Lammps_dumpreader(args.lammpstrj,m_array=1.) # m_array=mass1
    
    #移动盒子
    if args.box_mid:
        r.box_mid()

    ######################frame 0 read of skip
    if frame_strat == 0:
        condition = True 
    else:
        condition = False

    while r._read_single_frame(skip = not condition):
    
        if condition:
            #do analysis for every frame

            cm,rg = mf.get_mol_cm_rg(r._system,r.l_box,num_chain,num_beads)

            ## -rg
            if args.mol_rg:
                mol_rg.write(f'{r._timestep} ' + ' '.join(str(round(float(x),2)) for x in rg) +'\n')


            ##不同的找pair模式####################################
            if args.cluster_yzs:
                pair = mf.cluster_ysz(cm,num_chain,cut_off,r.l_box)
            if args.cluster_sy:
                pair,intrapair,outerpair = mf.cluster_sy(r._system,num_chain,num_beads,cut_off,r.l_box,r._type_index[1],r._type_index[2])
                bond_number.write(f'{r._timestep} {outerpair} {np.sum(intrapair)} {outerpair+np.sum(intrapair)} {len(r._type_index[1])+len(r._type_index[2])-2*(outerpair+np.sum(intrapair))}\n')
            #Union#######################################################
            cluster_list = mf.find_cluster(pair,num_chain)
            #aa = len(cluster_list)
            #print(aa)
            ###最大的cluster质心
            shift,_ = mf.cal_cm_rg(cm[cluster_list[0]][:,0],cm[cluster_list[0]][:,1],cm[cluster_list[0]][:,2],r.l_box,np.ones((len(cluster_list[0]),),dtype=int))
            
            for m in range(len(shift)):
                shift[m] = shift[m]-r.l_box[0]/2

            #shift[i for i in range(3)] = shift[i for i in range(3)]
            ###平移原子
            if args.cluster_mid:
                r.shift_pos(shift * (-1))
            # -cv
            if args.cluster_volum:
                cluster_size = collections.defaultdict(int)
                large_cluster_size = 0
                for i in cluster_list:
                    if len(i)>9:
                        large_cluster_size += len(i)
                        cluster_size[10] += 1
                    else:
                        cluster_size[len(i)] += 1
                cluster_volum.write(f'{r._timestep} {cluster_size[1]} {cluster_size[2]} {cluster_size[3]} {cluster_size[4]} {cluster_size[5]} {cluster_size[6]} {cluster_size[7]} {cluster_size[8]} {cluster_size[9]} {cluster_size[10]} {large_cluster_size}\n')
            ##  cluster -c
            if args.cluster:
                cluster.write(f'timestep = {r._timestep}\n')      
                averg = 0.0
                drop_num = 0
                drop_len = 0
            #    print (cluster_list)
                for m in cluster_list:
                    if len(m) > 1:
                        outter = 0
                        mc_number = len(m)
                        ty1 = [7,15,23]
                        ty2 = [31,39,47]
                        for i in range(mc_number-1):
                            ii = m[i]
                            for j in range(i+1,mc_number):
                                jj = m[j]
                                if (ii,jj) in pair: 
                                    outter += pair[(ii,jj)]
                        # tem_list = mf.cluster_to_ids(i,num_beads)
                #       _,tem_rg = mf.cal_cm_rg(r._system['x'][tem_list],r._system['y'][tem_list],r._system['z'][tem_list],r.l_box,r._system['mass'][tem_list])
                        cluster.write(f'{outter} {len(m)} ( ' + ' '.join(str(x) for x in m)+' )\n')
                    #    averg += tem_rg
                    #    drop_num += 1
                    #    drop_len += len(i)
                #if drop_num:
                #    drop.write(f'{r._timestep} {averg/drop_num} {drop_num} {drop_len}\n')
                #else:
                #    drop.write(f'{r._timestep} {drop_num} {drop_num} {drop_len}\n')
            # max_cluster
            if args.max_cluster:
                max_cluster.write(f'{r._timestep} {len(cluster_list[0])} ( ' + ' '.join(str(x) for x in cluster_list[0])+' )\n')

            if args.max_interbonds:
                inner = np.sum(intrapair[cluster_list[0]])
                outter = 0
                aa = 0
                internumber1 = 0
                internumber2 = 0
                internumber3 = 0
                internumber4 = 0
                internumber5 = 0
                internumber6 = 0
                internumber7 = 0
                internumber8 = 0
                internumber9 = 0
                internumber10 = 0
                internumber11 = 0
                internumber12 = 0
                internumber13 = 0
                internumber14 = 0
                internumber15 = 0
                internumber16 = 0
                # cdx = 0
                # cdy = 0
                # cdz = 0
                mc_number = len(cluster_list[0])
                ty1 = [7,15,23]
                ty2 = [31,39,47]
                for i in range(mc_number-1):
                    ii = cluster_list[0][i]
                    for j in range(i+1,mc_number):
                        jj = cluster_list[0][j]

                        if (ii,jj) in pair:
                            for k in range(7,48,8):
                                for kk in range(7,48,8):
                                    if (k in ty1 and kk in ty2) or (k in ty2 and kk in ty1):
                                        max_dx = abs(r._system['x'][((ii)*50+k)]-r._system['x'][((jj)*50+kk)])
                                        max_dy = abs(r._system['y'][((ii)*50+k)]-r._system['y'][((jj)*50+kk)])
                                        max_dz = abs(r._system['z'][((ii)*50+k)]-r._system['z'][((jj)*50+kk)])
                                        distance = max_dx**2 + max_dy**2 + max_dz**2
                                     #   print((ii)*50+k,(jj)*50+kk,r._system['z'][((ii)*50+k)],r._system['z'][((jj)*50+k)],np.sqrt(distance))
                                        if distance <= 3.5*3.5:
                                            cdx = abs(r._system['x'][((ii)*50+k)]-225)
                                            cdy = abs(r._system['y'][((ii)*50+k)]-225)
                                            cdz = abs(r._system['z'][((ii)*50+k)]-225)
                                            c_distance = np.sqrt(cdx**2 + cdy**2 + cdz**2)
                                            if c_distance <= 5:
                                                internumber1 += 1
                                            elif c_distance > 5 and c_distance <= 10:
                                                internumber2 += 1
                                            elif c_distance > 10 and c_distance <= 15:
                                                internumber3 += 1
                                            elif c_distance > 15 and c_distance <= 20:
                                                internumber4 += 1
                                            elif c_distance > 20 and c_distance <= 25:
                                                internumber5 += 1
                                            elif c_distance > 25 and c_distance <= 30:
                                                internumber6 += 1
                                            elif c_distance > 30 and c_distance <= 35:
                                                internumber7 += 1
                                            elif c_distance > 35 and c_distance <= 40:
                                                internumber8 += 1
                                            elif c_distance > 40 and c_distance <= 45:
                                                internumber9 += 1
                                            elif c_distance > 45 and c_distance <= 50:
                                                internumber10 += 1   
                                            elif c_distance > 50 and c_distance <= 55:
                                                internumber11 += 1 
                                            elif c_distance > 55 and c_distance <= 60:
                                                internumber12 += 1
                                            elif c_distance > 60 and c_distance <= 65:
                                                internumber13 += 1
                                            elif c_distance > 65 and c_distance <= 70:
                                                internumber14 += 1 
                                            elif c_distance > 70 and c_distance <= 75:
                                                internumber15 += 1 
                                            elif c_distance > 75 and c_distance <= 80:
                                                internumber16 += 1 
                            outter += pair[(ii,jj)]
            #    max_interbonds.write(f'{r._timestep} {outter} {inner} {outter+inner} {2*(len(r._type_index[1])//num_beads)*mc_number - 2*outter - 2*inner}\n')
                max_interbonds.write(f'{r._timestep} {outter} {internumber1/523.3} {internumber2/3663.3} {internumber3/9943.3} {internumber4/19363.3} {internumber5/31923.3} {internumber6/47623.3} {internumber7/66463.3} {internumber8/88443.3} {internumber9/113563.3} {internumber10/141823.3} {internumber11/173223.3} {internumber12/207763.3} {internumber13/245443.3} {internumber14/286263.3} {internumber15/3302223.3} {internumber16/377323.3}\n')
            if args.output:
                output.write(r.out_put_lammpstrj())
        # reader condition
        
        condition = True  if (not args.read) or ((frame_strat<= r.frame_num <=frame_end) and ((r.frame_num - frame_strat) % frame_step == 0)) else False

    del r


    
    end=datetime.datetime.now()
    print(f'start at {start}')
    print(f'end   at {end}')
    print(f'total wall time is  {end-start}')

    
if __name__ == "__main__":

    #open file
    # -o
    if args.output:
        output = open(args.output,'w+')
    # -rg
    if args.mol_rg:
        mol_rg = open(args.lammpstrj + '_rg.txt','w+')
    # -c
    if args.cluster:
        cluster = open('cluster.txt','w+')
        drop = open('drop.txt','w+')
    # -cv
    if args.cluster_volum:
        cluster_volum = open('cluster_volum.txt','w+')
    # -mc
    if args.max_cluster:
        max_cluster = open('max_cluster.txt','w+')
    # -mcb
    if args.max_interbonds:
        max_interbonds = open('max_interbonds.txt','w+')
    # -sy
    if args.cluster_sy:
        bond_number = open('bond_number.txt','w+')
    
    ##############################
    main()
    ##############################

    #close file 
    if args.output:
        output.close()
    if args.mol_rg:
        mol_rg.close()
    if args.max_cluster:
        max_cluster.close()
    if args.cluster_sy:
        bond_number.close()
    if args.max_interbonds:
        max_interbonds.close()
    if args.cluster:
        cluster.close()
        drop.close()
    if args.cluster_volum:
        cluster_volum.close()


