

import numpy as np
import random
import matplotlib.pyplot as plt

# Initialize the values of parameters and coefficients.
# cycle time
# ct=
# number of tasks
# nt=
# maximum number of workstations
# nw=
# Define precedence relationship:
# p = [[0]*nt for i in range(nt)]
# if task i need task j, p[i][j]=1.


# cost of robots, cobots
# costr=
# costc=
# number of available cobots, robots
# qrr=
# qcc=



# processing time of baseline with only human operations
# proc=

# coefficient characterizing the processing time change with robot

# trh=
# tre=
# trd=
# trp=
# tri=

# coefficient characterizing the processing time change with cobot or robot: rhocc, rhorr
# rhocc=
# rhorr=

# coefficient characterizing the strain factor change with cobot or robot:
# betapp=
# betaii=



# RSI multiplier M_I of tasks

def r2im(r):
    if r<= 0.4:
        return (30*r**3-15.6*r**2+13*r+0.4)
    else:
        return (36*r**3-33.3*r**2+24.77*r-1.86)

# RSI multiplier M_P of tasks

def r2pm(r):
    if r <= 0:
        return (1.2*np.exp(-0.009*r)-0.2)
    elif r<= 30:
        return (1)
    else:
        return (1+0.00028*(r-30)**2)

# RSI multiplier M_D of tasks

def r2dm(r):
    if r <= 60:
        return (0.45+0.31*r)
    else:
        return (19.17*np.log(r)-59.44)

# RSI multiplier M_E of tasks

def r2em(r):
    if r <= 90:
        return (0.1+0.25*r)
    else:
        return (0.00334*r**1.96)


# RSI multiplier M_H of tasks

def r2hm(r):
    if r<= 0.05:
        return (0.2)
    else:
        return (0.042*r+0.090*np.log(r)+0.477)

# Define stain constraint
# K=

# Revised strain index (RSI)
def rsi(e,h,d,p,i):
    return(r2im(i)*r2em(e)*r2dm(d)*r2pm(p)*r2hm(h))

# Composite strain index (COSI)
def cosi(rsi_list,e_list):
    e_list.sort(key=dict(zip(e_list,rsi_list)).get,reverse=True)
    rsi_list.sort(reverse=True)
    firsi_list=[rsi_list[i]/e_list[i] for i in range(len(rsi_list))]
    acc_e=[sum(e_list[:(i+1)]) for i in range(len(rsi_list))]
    delta_em=[r2em(acc_e[i+1])-r2em(acc_e[i]) for i in range(len(rsi_list)-1)]
    delta_rsi=[firsi_list[i+1]*delta_em[i] for i in range(len(rsi_list)-1)]
    cosi=rsi_list[0]+sum(delta_rsi)
    return (cosi)

# Cumulative strain index (CUSI)
def cusi(cosi_list,h_list):
    h_list.sort(key=dict(zip(h_list,cosi_list)).get,reverse=True)
    cosi_list.sort(reverse=True)
    hicosi_list=[cosi_list[i]/h_list[i] for i in range(len(cosi_list))]
    acc_h=[sum(h_list[:(i+1)]) for i in range(len(cosi_list))]
    delta_hm=[r2hm(acc_h[i+1])-r2em(acc_h[i]) for i in range(len(cosi_list)-1)]
    delta_cosi=[hicosi_list[i+1]*delta_hm[i] for i in range(len(cosi_list)-1)]
    cusi=cosi_list[0]+sum(delta_cosi)
    return (cusi)

# random generate solution
def generate_solution(qr,qc,nt):
    seq=random.sample(range(0,nt),nt)
    alt=['h']*nt+['r']*qrr+['c']*qcc
    alt=random.sample(alt,nt)
    return [seq,alt]
dic_rho={'h':rhohh,'c':rhocc,'r':rhorr}

# transform solution representation into workstation representation
def split_workstation(sol):
    seq, alt=sol
    work=[1]*nt
    t=0
    t_list=[]
    w=0
    new_proc=[proc[i]*dic_rho[alt[i]][i] for i in range(nt)]
    for i in range(nt):
        t=t+new_proc[seq[i]]
        if t<=ct:
            work[i]=w
        else:
            w=w+1
            t_list.append(t-new_proc[seq[i]])
            t=new_proc[seq[i]]
            work[i]=w
    t_w=[work[seq[i]] for i in range(nt)]
    t_list.append(t)
    return [t_w,t_list]

# check if meet workstation number requirement
def check_workstation(t_w):
    if max(t_w)+1<=nw:
        return True
    else:
        return False

# check if meet precedence requirement
def check_precedence(sol):
    seq, alt=sol
    ord=np.argsort(seq)
    pp=True
    for i in range(nt):
        for j in range(nt):
            if (p[i][j]==1) & (ord[i]<ord[j]):
                pp=False
    return pp

# check if cosi meet requirement
def check_cosi(sol):
    seq, alt=sol
    t_w,t_list=split_workstation(sol)
    all_cosi=[]
    ss=True
    for i in range(max(t_w)+1):
        e_list=[]
        h_list=[]
        i_list=[]
        p_list=[]
        d_list=[]
        rsi_list=[]
        for j in range(nt):
            if (t_w[j]==i) & (alt[j]!='r'):
                e_list.append(tre[alt[j]][j])
                h_list.append(trh[alt[j]][j])
                d_list.append(trd[alt[j]][j])
                i_list.append(tri[alt[j]][j])
                p_list.append(trp[alt[j]][j])
                rsi_list.append(rsi(tre[alt[j]][j],trh[alt[j]][j],trd[alt[j]][j],trp[alt[j]][j],tri[alt[j]][j]))
        if len(rsi_list)>0:
            cosi_i=cosi(rsi_list,e_list)
        else:
            cosi_i=0
        if cosi_i>K:
            ss=False
        all_cosi.append(cosi_i)
    return [all_cosi,ss]

# crossover operator 1
def crossover1(sol,sol1):
    seq1,alt1=sol1
    seq,alt=sol
    new_seq=seq[:]
    new_alt=alt[:]
    positions=random.sample(range(nt),2)
    positions.sort()
    pos1,pos2=positions
    new_seq[pos1:pos2]=[posi for posi in seq1 if posi in seq[pos1:pos2]]
    new_alt[pos1:pos2]=[alt1[i] for i in range(nt) if seq1[i] in seq[pos1:pos2]]
    return [new_seq,new_alt]
# crossover operator 2
def crossover2(sol,sol1):
    seq1,alt1=sol1
    seq,alt=sol
    new_seq=seq[:]
    new_alt=alt[:]
    pos=random.sample(range(nt),1)[0]
    new_seq[pos:]=[posi for posi in seq1 if posi in seq[pos:]]
    new_alt[pos:]=[alt1[i] for i in range(nt) if seq1[i] in seq[pos:]]
    return [new_seq,new_alt]
# crossover operator 3
def crossover3(sol,sol1):
    seq1,alt1=sol1
    seq,alt=sol
    new_seq=seq[:]
    new_alt=alt[:]
    positions=random.sample(range(nt),2)
    positions.sort()
    pos1,pos2=positions
    new_seq[pos2:]=[posi for posi in seq1 if posi in seq[pos2:]]
    new_seq[:pos1]=[posi for posi in seq1 if posi in seq[:pos1]]
    new_alt[pos2:]=[alt1[i] for i in range(nt) if seq1[i] in seq[pos2:]]
    new_alt[:pos1]=[alt1[i] for i in range(nt) if seq1[i] in seq[:pos1]]
    return [new_seq,new_alt]
# crossover operator 4
def crossover4(sol,sol1):
    seq1,alt1=sol1
    seq,alt=sol
    new_seq=seq[:]
    new_alt=alt[:]
    nr=random.sample(range(nt+1),1)[0]
    positions=random.sample(range(nt),nr)
    m_seq=[]
    j=0
    for i in range(nt):
        if i not in positions:
            m_seq.append(seq[i])
    m1_seq=[seq1[i] for i in range(nt) if seq1[i] in m_seq]
    m1_alt=[alt1[i] for i in range(nt) if seq1[i] in m_seq]
    for i in range(nt):
        if i not in positions:
            new_seq[i]=m1_seq[j]
            new_alt[i]=m1_alt[j]
            j=j+1
    return [new_seq,new_alt]
# crossover operator 5
def crossover5(sol,sol1):
    seq1,alt1=sol1
    seq,alt=sol
    new_seq=seq[:]
    new_alt=alt[:]
    pos=random.sample(range(nt),1)[0]
    pos1=pos
    m_seq=[]
    m_alt=[]
    j=0
    temp=seq1[pos]
    temp_alt=alt1[pos]
    pos=seq.index(temp)
    new_seq[pos]=temp
    new_alt[pos]=temp_alt
    while pos1!=pos:
        temp=seq1[pos]
        temp_alt=alt1[pos]
        pos=seq.index(temp)
        new_seq[pos]=temp
        new_alt[pos]=temp_alt
    return [new_seq,new_alt]

# mutation operator 1
def mutation1(sol):
    seq,alt=sol
    positions=random.sample(range(nt),2)
    positions.sort()
    pos1,pos2=positions
    new_seq=seq[:]
    new_seq[pos1]=seq[pos2]
    new_seq[pos2]=seq[pos1]
    new_alt=alt[:]
    new_alt[pos1]=alt[pos2]
    new_alt[pos2]=alt[pos1]
    return [new_seq,new_alt]

# mutation operator 2
def mutation2(sol):
    seq,alt=sol
    positions=random.sample(range(nt),2)
    positions.sort()
    pos1,pos2=positions
    new_seq=seq[:]
    new_seq[pos1:pos2]=seq[(pos1+1):(pos2+1)]
    new_seq[pos2]=seq[pos1]
    new_alt=alt[:]
    new_alt[pos1:pos2]=alt[(pos1+1):(pos2+1)]
    new_alt[pos2]=alt[pos1]
    return [new_seq,new_alt]

# mutation operator 3
def mutation3(sol):
    seq,alt=sol
    positions=random.sample(range(nt),2)
    positions.sort()
    pos1,pos2=positions
    new_seq=seq[:]
    new_seq[(pos1+1):(pos2+1)]=seq[pos1:pos2]
    new_seq[pos1]=seq[pos2]
    new_alt=alt[:]
    new_alt[(pos1+1):(pos2+1)]=alt[pos1:pos2]
    new_alt[pos1]=alt[pos2]
    return [new_seq,new_alt]

# switch operator
def switch_new(sol,from_alt='h',to_alt='r'):
    seq,alt=sol
    from_indices=[i for i in range(nt) if alt[i]==from_alt]
    if len(from_indices)==0:
        return sol
    pos1=random.sample(from_indices,1)[0]
    new_alt=alt[:]
    new_alt[pos1]=to_alt
    return [seq,new_alt]
# define coefficients in reward function
# eta_ct=
# eta_w=
# pen_st=
# pen_prec=
# pen_c=
# pen_r=
# reward calculation
def get_reward(sol):
    reward = 0
    t_w, ct_list = split_workstation(sol)
    num_workstation = max(t_w) + 1
    actual_ct = max(ct_list)
    all_cosi, ss = check_cosi(sol)
    pp = check_precedence(sol)
    if ss == False: reward = reward - pen_st
    if pp == False: reward = reward - pen_prec
    alt = sol[1]
    qc = np.sum([alt_i == 'c' for alt_i in alt])
    qr = np.sum([alt_i == 'r' for alt_i in alt])
    reward = reward - qc * costc - qr * costr
    if qc > qcc:
        reward = reward - pen_c
    if qr > qrr:
        reward = reward - pen_r
    reward = reward - actual_ct * eta_ct - sum(all_cosi) - num_workstation * eta_w
    return reward

# check the feasibility of the solution

def if_valid(sol):
    vv = True
    all_cosi, ss = check_cosi(sol)
    pp = check_precedence(sol)
    alt = sol[1]
    if np.sum([alt_i == 'c' for alt_i in alt]) > qcc:
        vv = False
    if np.sum([alt_i == 'r' for alt_i in alt]) > qrr:
        vv = False
    if ss == False: vv = False
    if pp == False: vv = False
    return vv

# check if the iteration generate useful solution

def if_terminate(sol_set, new_sol_set):
    exist_reward = [get_reward(sol_i) for sol_i in sol_set]
    new_reward = [get_reward(sol_i) for sol_i in new_sol_set]
    if max(new_reward) <= min(exist_reward):
        return True
    else:
        return False

# check if the iteration improve the solution set

def if_improve(sol_set, new_sol_set):
    exist_reward = [get_reward(sol_i) for sol_i in sol_set]
    new_reward = [get_reward(sol_i) for sol_i in new_sol_set]
    ii = 1
    for i in range(len(exist_reward)):
        if new_reward[i] > exist_reward[i]:
            ii = 0
    return ii

# function to update solution set

def update_sol_set(sol_set, new_sol_set, ps):
    all_sol = sol_set[:]
    for sol_i in new_sol_set:
        if sol_i not in all_sol:
            all_sol.append(sol_i)
    all_reward = np.array([get_reward(sol_i) for sol_i in all_sol])
    order = all_reward.argsort()
    ranked_sol = [all_sol[i] for i in order]
    ranked_sol.reverse()
    return ranked_sol[:ps]

def argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

# define simulation parameters
# num_sim=
# max_run=

# define size of solution set
# ps=

# initialize solution set
initial_set=[]
while (len(initial_set)<=ps):
    rand_sol=generate_solution(qrr,qcc,nt)
    if (rand_sol not in initial_set):
        initial_set.append(rand_sol)

# define dynamic hyperparameter
def dynamic_epsilon(iter):
    return 0.2 - iter / max_run * 0.2
def dynamic_CR(iter):
    return 0.7 - 0.35 * iter/ max_run
# define Q-learning hyperparameter
# discount_factor =
# learning_rate =

# run Q_EA
obj_list_list = []
for nr in range(num_sim):
    state_list = []
    action_list = []
    sol_set = initial_set

    # states
    states = ['nonimprove', 'improve']
    current_state = 0
    state_value = {'nonimprove': 0, 'improve': 1}
    q_cr = np.zeros((2, 5))
    q_mu = np.zeros((2, 3))
    q_sw = np.zeros((2, 6))
    SR = 0.3
    num_run = 0
    old_sol_set = sol_set[:]
    visited_sol = []
    obj_list = []
    while (num_run < max_run):
        epsilon = dynamic_epsilon(num_run)
        CR = dynamic_CR(num_run)
        MR = 0.7 - CR
        old_sol_set = sol_set[:]
        num_run = num_run + 1
        # cross over phase
        temp_set = []
        current_reward_bound = np.mean([get_reward(sol_i) for sol_i in sol_set])
        obj_list.append(max([get_reward(sol_i) for sol_i in sol_set]))
        if random.uniform(0, 1) < CR:
            ana_sol = random.choice(sol_set)
            if random.uniform(0, 1) < epsilon:
                action_index = argmax(q_cr[current_state, :])
            else:
                action_index = random.sample(range(5), 1)[0]
            for sol_i in sol_set:
                if action_index == 0:
                    new_cr_sol = crossover1(sol_i, ana_sol)
                elif action_index == 1:
                    new_cr_sol = crossover2(sol_i, ana_sol)
                elif action_index == 2:
                    new_cr_sol = crossover3(sol_i, ana_sol)
                elif action_index == 3:
                    new_cr_sol = crossover4(sol_i, ana_sol)
                elif action_index == 4:
                    new_cr_sol = crossover5(sol_i, ana_sol)
                temp_set.append(new_cr_sol)
                if (new_cr_sol not in visited_sol):
                    visited_sol.append(new_cr_sol)
            old_state = current_state
            current_state = if_improve(sol_set, temp_set)
            state_list.append(current_state)
            action_list.append([1, action_index])
            old_q_cr = q_cr[old_state, action_index]
            temp_diff = np.mean(
                [get_reward(sol_i) for sol_i in temp_set]) - current_reward_bound + discount_factor * np.max(
                q_cr[current_state]) - old_q_cr
            q_cr[old_state, action_index] = old_q_cr + learning_rate * temp_diff
            sol_set = update_sol_set(sol_set, temp_set, ps)

        # mutation phase
        temp_set = []
        if random.uniform(0, 1) < MR:
            if random.uniform(0, 1) < epsilon:
                action_index = argmax(q_mu[current_state, :])
            else:
                action_index = random.sample(range(3), 1)[0]
            for sol_i in sol_set:
                if action_index == 0:
                    new_mu_sol = mutation1(sol_i)
                elif action_index == 1:
                    new_mu_sol = mutation2(sol_i)
                elif action_index == 2:
                    new_mu_sol = mutation3(sol_i)
                if (new_mu_sol not in visited_sol):
                    visited_sol.append(new_mu_sol)
                temp_set.append(new_mu_sol)
            old_state = current_state
            current_state = if_improve(sol_set, temp_set)
            action_list.append([2, action_index])
            state_list.append(current_state)
            old_q_mu = q_mu[old_state, action_index]
            temp_diff = np.mean(
                [get_reward(sol_i) for sol_i in temp_set]) - current_reward_bound + discount_factor * np.max(
                q_mu[current_state]) - old_q_mu
            q_mu[old_state, action_index] = old_q_mu + learning_rate * temp_diff
            sol_set = update_sol_set(sol_set, temp_set, ps)

        # switch phase
        temp_set = []
        if random.uniform(0, 1) < SR:
            if random.uniform(0, 1) < epsilon:
                action_index = argmax(q_mu[current_state, :])
            else:
                action_index = random.sample(range(6), 1)[0]
            for sol_i in sol_set:
                if action_index == 0:
                    new_sw_sol = switch_new(sol_i, 'h', 'r')
                elif action_index == 1:
                    new_sw_sol = switch_new(sol_i, 'h', 'c')
                elif action_index == 2:
                    new_sw_sol = switch_new(sol_i, 'r', 'c')
                elif action_index == 3:
                    new_sw_sol = switch_new(sol_i, 'r', 'h')
                elif action_index == 4:
                    new_sw_sol = switch_new(sol_i, 'c', 'r')
                elif action_index == 5:
                    new_sw_sol = switch_new(sol_i, 'c', 'h')
                if (new_sw_sol not in visited_sol):
                    visited_sol.append(new_sw_sol)
                temp_set.append(new_sw_sol)
            old_state = current_state
            current_state = if_improve(sol_set, temp_set)
            action_list.append([3, action_index])
            state_list.append(current_state)
            old_q_sw = q_sw[old_state, action_index]
            temp_diff = np.mean(
                [get_reward(sol_i) for sol_i in temp_set]) - current_reward_bound + discount_factor * np.max(
                q_sw[current_state]) - old_q_sw
            q_sw[old_state, action_index] = old_q_sw + learning_rate * temp_diff
            sol_set = update_sol_set(sol_set, temp_set, ps)
    print([get_reward(sol_i) for sol_i in sol_set])
    print([if_valid(sol_i) for sol_i in sol_set])
    print([max(split_workstation(sol_i)[0]) + 1 for sol_i in sol_set])
    print(len(visited_sol))
    print(sol_set)
    obj_list_list.append(obj_list)