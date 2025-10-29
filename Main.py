# ====== imports ======
import math
from pyomo import environ as pyo
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ====== input data ======
units =  ['compressor house','furnace','reactor','separators','distil','reboiler',
          'methanol tank','h2 tank','co2 tank','battery','ctrlroom']

alpha = {'compressor house':30,'furnace':10,'reactor':20,'separators':10,'distil':15,'reboiler':10,
         'methanol tank':20,'h2 tank':75,'co2 tank':20,'battery':35,'ctrlroom':20}
beta  = {'compressor house':30,'furnace':10,'reactor':20,'separators':10,'distil':15,'reboiler':10,
         'methanol tank':30,'h2 tank':75,'co2 tank':20,'battery':25,'ctrlroom':20}

# ====== safety distances  ======
Demin = {}
def setD(i,j,val):
    Demin[i,j]=val; Demin[j,i]=val
for i in units:
    for j in units:
        if i!=j: Demin[i,j]=0.0

# 0: compressor house
setD('compressor house','furnace',150) #meters (m in safety radii)
setD('compressor house','reactor',150)
setD('compressor house','separators',150)
setD('compressor house','distil',150)
setD('compressor house','reboiler',150)
setD('compressor house','methanol tank',150)
setD('compressor house','h2 tank',150)
setD('compressor house','co2 tank',150)
setD('compressor house','battery',150)
setD('compressor house','ctrlroom',150)
# 1: furnace
setD('furnace','reactor',100)
setD('furnace','separators',100)
setD('furnace','distil',100)
setD('furnace','reboiler',100)
setD('furnace','methanol tank',100)
setD('furnace','h2 tank',100)
setD('furnace','co2 tank',100)
setD('furnace','battery',100)
setD('furnace','ctrlroom',100)

# 2: reactor
setD('reactor','furnace',130)
setD('reactor','separators',130)
setD('reactor','distil',130)
setD('reactor','reboiler',130)
setD('reactor','methanol tank',130)
setD('reactor','h2 tank',130)
setD('reactor','co2 tank',130)
setD('reactor','battery',130)
setD('reactor','ctrlroom',130)
# 3: separators
setD('separators','distil',130)
setD('separators','reboiler',130)
setD('separators','methanol tank',130)
setD('separators','h2 tank',130)
setD('separators','co2 tank',130)
setD('separators','battery',130)
setD('separators','ctrlroom',130)
# 4: distil
setD('distil','reboiler',120)
setD('distil','methanol tank',120)
setD('distil','h2 tank',120)
setD('distil','co2 tank',120)
setD('distil','battery',120)
setD('distil','ctrlroom',120)
# 5: reboiler
setD('reboiler','methanol tank',125)
setD('reboiler','h2 tank',125)
setD('reboiler','co2 tank',125)
setD('reboiler','battery',125)
setD('reboiler','ctrlroom',125)
# 6: methanol tank
setD('methanol tank','h2 tank',100.0)
setD('methanol tank','co2 tank',100.0)
setD('methanol tank','battery',100.0)
setD('methanol tank','ctrlroom',100.0)
# 7: h2 tank
setD('h2 tank','co2 tank',100)
setD('h2 tank','battery',100)
setD('h2 tank','ctrlroom',100)
# 8: co2 tank
setD('co2 tank','battery',90)
setD('co2 tank','ctrlroom',90)
# 9: battery
setD('battery','ctrlroom',80)

# ====== piping properties  ======
velocity = {(i,j):0.0 for i in units for j in units}
Q        = {(i,j):0.0 for i in units for j in units}   # volumetric m^3/s
npp      = {(i,j):0   for i in units for j in units}   # number of parallel pipes

def set_pair(mat,i,j,val): 
    mat[i,j]=val; mat[j,i]=val

set_pair(velocity, 'compressor house','furnace', 30.0)
set_pair(velocity, 'furnace','reactor', 30.0)
set_pair(velocity, 'reactor','separators', 30.0)
set_pair(velocity, 'separators','distil', 3.0)
set_pair(velocity, 'distil','reboiler', 3.0)
set_pair(velocity, 'reboiler','battery', 3.0)
set_pair(velocity, 'distil','methanol tank', 3.0)
set_pair(velocity, 'separators','compressor house', 30.0)
set_pair(velocity, 'h2 tank','compressor house', 30.0)
set_pair(velocity, 'co2 tank','compressor house', 30.0)

# --- volumetric flow Q (m^3/s) ---
set_pair(Q, 'compressor house', 'furnace',       116658.3/3600)
set_pair(Q, 'furnace','reactor',       116658.3/3600)
set_pair(Q, 'reactor','separators',     95448.9/3600)
set_pair(Q, 'separators','distil',         20143.6/3600)
set_pair(Q, 'distil','reboiler',        7409.2/3600)
set_pair(Q, 'reboiler','battery',        7409.2/3600)
set_pair(Q, 'distil','methanol tank',  12499.7/3600)
set_pair(Q, 'separators','compressor house', 90676.3/3600)
set_pair(Q, 'h2 tank','compressor house', 22381.6/3600)
set_pair(Q, 'co2 tank','compressor house',  3075.6/3600)

# --- number of parallel pipes ---
for a,b in [
    ('compressor house','furnace'),
    ('furnace','reactor'),
    ('reactor','separators'),
    ('separators','distil'),
    ('distil','reboiler'),
    ('reboiler','battery'),
    ('distil','methanol tank'),
    ('separators','compressor house'),
    ('h2 tank','compressor house'),
    ('co2 tank','compressor house'),
]:
    npp[a,b] = 1
    npp[b,a] = 1
# ====== cost constants ======
LAND_PRICE = 61.78   
BB          = 880
n_2         = 0.74
CEPCI_2006  = 499.6
CEPCI_2021  = 677.7
FX_rate     = 0.75
F_install   = 1.5

# ====== build pipe cost coefficients ======
c_per_m = {(i,j):0.0 for i in units for j in units}
pipe_pairs_raw = []
for j in units:
    for i in units:
        if units.index(j) > units.index(i) and npp[i,j]>0 and velocity[i,j]>0 and Q[i,j]>0:
            Qvol = Q[i,j]
            DIA = math.sqrt((4*Qvol) / (velocity[i,j]*math.pi))
            PC_bare = BB * (DIA**n_2) * (CEPCI_2021/CEPCI_2006) * FX_rate
            PC_installed = PC_bare * (1 + F_install)
            c_per_m[i,j] = PC_installed * npp[i,j]
            c_per_m[j,i] = c_per_m[i,j]
            pipe_pairs_raw.append((i,j))

# ====== pair sets ======
all_pairs  = [(ui, uj) for a, ui in enumerate(units) for uj in units[a+1:]] 
pipe_pairs = list(pipe_pairs_raw) 

# ====== model ======
m = pyo.ConcreteModel()
m.U      = pyo.Set(initialize=units)
m.Pall   = pyo.Set(initialize=all_pairs,  dimen=2)   
m.Ppipes = pyo.Set(initialize=pipe_pairs, dimen=2)  

# Vars
m.O  = pyo.Var(m.U, within=pyo.Binary)
m.x  = pyo.Var(m.U, within=pyo.NonNegativeReals)
m.y  = pyo.Var(m.U, within=pyo.NonNegativeReals)
m.l  = pyo.Var(m.U, within=pyo.NonNegativeReals)
m.d  = pyo.Var(m.U, within=pyo.NonNegativeReals)

# L1 distance helpers 
m.dx = pyo.Var(m.Pall, within=pyo.NonNegativeReals)
m.dy = pyo.Var(m.Pall, within=pyo.NonNegativeReals)
m.t  = pyo.Var(m.Pall, within=pyo.NonNegativeReals)

# ---- Non-overlap   ----
Mbig = 10000
m.bx = pyo.Var(m.Pall, within=pyo.Binary)  
m.sx = pyo.Var(m.Pall, within=pyo.Binary)  
m.sy = pyo.Var(m.Pall, within=pyo.Binary) 

# Land bounding box 
m.Xmin = pyo.Var(within=pyo.NonNegativeReals)
m.Ymin = pyo.Var(within=pyo.NonNegativeReals)
m.Xmax = pyo.Var(within=pyo.NonNegativeReals)
m.Ymax = pyo.Var(within=pyo.NonNegativeReals)
m.Xmin.fix(0.0)
m.Ymin.fix(0.0)

# Orientation & geometry
m.c_orient_l = pyo.Constraint(m.U, rule=lambda m,i: m.l[i] == alpha[i]*m.O[i] + beta[i]*(1 - m.O[i]))
m.c_orient_d = pyo.Constraint(m.U, rule=lambda m,i: m.d[i] == alpha[i] + beta[i] - m.l[i])

# Keep inside land box
m.c_box1 = pyo.Constraint(m.U, rule=lambda m,i: m.x[i] - 0.5*m.l[i] >= m.Xmin)
m.c_box2 = pyo.Constraint(m.U, rule=lambda m,i: m.y[i] - 0.5*m.d[i] >= m.Ymin)
m.c_box3 = pyo.Constraint(m.U, rule=lambda m,i: m.x[i] + 0.5*m.l[i] <= m.Xmax)
m.c_box4 = pyo.Constraint(m.U, rule=lambda m,i: m.y[i] + 0.5*m.d[i] <= m.Ymax)

# |x_i - x_j|, |y_i - y_j| (only for t)
m.c_dx1 = pyo.Constraint(m.Pall, rule=lambda m,i,j: m.dx[i,j] >=  m.x[i] - m.x[j])
m.c_dx2 = pyo.Constraint(m.Pall, rule=lambda m,i,j: m.dx[i,j] >=  m.x[j] - m.x[i])
m.c_dy1 = pyo.Constraint(m.Pall, rule=lambda m,i,j: m.dy[i,j] >=  m.y[i] - m.y[j])
m.c_dy2 = pyo.Constraint(m.Pall, rule=lambda m,i,j: m.dy[i,j] >=  m.y[j] - m.y[i])
m.c_t   = pyo.Constraint(m.Pall, rule=lambda m,i,j: m.t[i,j]  ==  m.dx[i,j] + m.dy[i,j])

# --- Required center offsets to ensure edge-to-edge + safety ---
def req_x(m,i,j):  
    return 0.5*(m.l[i] + m.l[j]) + Demin.get((i,j), 0.0)

def req_y(m,i,j):  
    return 0.5*(m.d[i] + m.d[j]) + Demin.get((i,j), 0.0)


m.sep_x1 = pyo.Constraint(m.Pall, rule=lambda m,i,j:
    m.x[i] - m.x[j] >= req_x(m,i,j) - Mbig*(1 - m.bx[i,j]) - Mbig*(1 - m.sx[i,j]))
m.sep_x2 = pyo.Constraint(m.Pall, rule=lambda m,i,j:
    m.x[j] - m.x[i] >= req_x(m,i,j) - Mbig*(1 - m.bx[i,j]) - Mbig*(    m.sx[i,j]))

m.sep_y1 = pyo.Constraint(m.Pall, rule=lambda m,i,j:
    m.y[i] - m.y[j] >= req_y(m,i,j) - Mbig*(    m.bx[i,j]) - Mbig*(1 - m.sy[i,j]))
m.sep_y2 = pyo.Constraint(m.Pall, rule=lambda m,i,j:
    m.y[j] - m.y[i] >= req_y(m,i,j) - Mbig*(    m.bx[i,j]) - Mbig*(    m.sy[i,j]))

# Fixing a unit
m.x['methanol tank'].fix(alpha['methanol tank']/2)
m.y['methanol tank'].fix(beta['methanol tank']/2)

# ----- Land CAPEX: piecewise-linear outer-approx of area Z >= U*V -----
m.Uw = pyo.Var(within=pyo.NonNegativeReals)  
m.Vh = pyo.Var(within=pyo.NonNegativeReals)  
m.Za = pyo.Var(within=pyo.NonNegativeReals)  

m.c_Uw = pyo.Constraint(expr=m.Uw == m.Xmax - m.Xmin)
m.c_Vh = pyo.Constraint(expr=m.Vh == m.Ymax - m.Ymin)

U_ub = sum(beta[u] for u in units) + 500.0
V_ub = sum(beta[u] for u in units) + 500.0

u_pts = [0, 100, 300, 600, 900, 1200, min(900.0, U_ub), U_ub]
v_pts = [0, 100, 300, 600, 900, 1200, min(900.0, V_ub), V_ub]
m.area_cuts = pyo.ConstraintList()
for uk in u_pts:
    for vk in v_pts:
        m.area_cuts.add( m.Za >= vk*m.Uw + uk*m.Vh - uk*vk )

m.land_capex = pyo.Expression(expr=LAND_PRICE * m.Za)

# ----- Pipe CAPEX -----
m.pipe_capex = pyo.Expression(expr=sum(c_per_m[i,j]*m.t[i,j] for (i,j) in m.Ppipes))

# Objective
m.obj = pyo.Objective(expr=m.land_capex + m.pipe_capex, sense=pyo.minimize)

# ====== solve  ======
solver = pyo.SolverFactory('gams:cplex')

solver.options['add_options'] = [
    'option reslim=300;',   # time limit (seconds)
    'option optcr=0.02;',   # relative optimality gap (2%)
]

res = solver.solve(m, tee=True)

# ====== report ======
print("Total CAPEX =", pyo.value(m.obj))
print("  Land CAPEX:", pyo.value(m.land_capex))
print("  Pipe CAPEX:", pyo.value(m.pipe_capex))
for u in units:
    print(f"{u:16s}  x={pyo.value(m.x[u]):8.2f}  y={pyo.value(m.y[u]):8.2f}  "
          f"l={pyo.value(m.l[u]):6.2f}  d={pyo.value(m.d[u]):6.2f}")
print("Bounding box:",
      "xmin", pyo.value(m.Xmin), "xmax", pyo.value(m.Xmax),
      "ymin", pyo.value(m.Ymin), "ymax", pyo.value(m.Ymax))


De = {
    'compressor house': 150,
    'furnace'         : 100,
    'reactor'         : 130,
    'separators'      : 130,
    'distil'          : 120,
    'reboiler'        : 125,
    'methanol tank'   : 100,
    'h2 tank'         : 100,
    'co2 tank'        : 90,
  
}

# ====== plot layout  ======
xmin = pyo.value(m.Xmin); xmax = pyo.value(m.Xmax)
ymin = pyo.value(m.Ymin); ymax = pyo.value(m.Ymax)
centers = {u: (pyo.value(m.x[u]), pyo.value(m.y[u])) for u in units}
sizes   = {u: (pyo.value(m.l[u]), pyo.value(m.d[u])) for u in units}

fig, ax = plt.subplots(figsize=(10, 8))
ax.set_aspect('equal')


bb_rect = mpatches.Rectangle((xmin, ymin), max(0.1, xmax - xmin), max(0.1, ymax - ymin),
                             fill=False, linewidth=2)
ax.add_patch(bb_rect)


for u in units:
    x, y = centers[u]; l, d = sizes[u]
    left = x - 0.5*l; bottom = y - 0.5*d


    rect = mpatches.Rectangle((left, bottom), l, d, fill=True, alpha=0.20, linewidth=2)
    ax.add_patch(rect)
    ax.annotate(u, (x, y), xytext=(3, 3), textcoords="offset points", fontsize=8)


    if u in De:
        r = De[u] + 0.5*max(l, d)  
        circ = mpatches.Circle((x, y), r, fill=False, linestyle='--', linewidth=1.5, alpha=0.6)
        ax.add_patch(circ)


for (i, j) in m.Ppipes:
    xi, yi = centers[i]; xj, yj = centers[j]
    ax.plot([xi, xj], [yi, yj])

pad_x = 0.1*(max(1.0, xmax - xmin))
pad_y = 0.1*(max(1.0, ymax - ymin))
ax.set_xlim((xmin - pad_x, xmax + pad_x))
ax.set_ylim((ymin - pad_y, ymax + pad_y))
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
ax.set_title(f"Plant layout — Total CAPEX = {pyo.value(m.obj):.2f}")


unit_patch   = mpatches.Patch(facecolor='C0', alpha=0.20, label='Unit footprint')
safety_patch = mpatches.Patch(fill=False, edgecolor='k', linestyle='--', label='Safety radius (from edge)')
ax.legend(handles=[unit_patch, safety_patch], loc='upper left')

plt.show()
