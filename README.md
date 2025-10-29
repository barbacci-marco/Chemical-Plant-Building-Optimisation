# Plant Layout Optimisation (Pyomo)

<figure>
  <img src="Figures/Final Optimisation Figure.png" alt="Plant layout with safety radii and pipes" width="700">
</figure>

## What this repo does

This script builds and solves a **mixed-integer nonlinear (MILP)** layout model for a process plant.
It chooses **x–y positions and orientations** of rectangular units to:

* minimise **total CAPEX = land cost + piping cost**, subject to
* **non-overlap** of units,
* **pairwise safety distances** between units (edge-to-edge),
* a variable **land bounding box** (the plot you must purchase).

It then **plots the layout**, including **unit footprints**, **optional safety radii (from unit edges)**, and **pipe polylines** between connected units.

---

## Key assumptions

* Units are **axis-aligned rectangles** with two possible orientations: `alpha[u] × beta[u]` or `beta[u] × alpha[u]`.
* Safety distance `Demin[i,j]` is an **extra margin** added **beyond rectangular edges** (implemented as a minimum **centre-to-centre** offset).
* Pipes are costed using **Manhattan (L1) distance** between unit centres. 
* The land is a **rectangle** `[Xmin, Xmax] × [Ymin, Ymax]`; here `Xmin=Ymin=0` are fixed to avoid drift. The **area** is modelled with a **piecewise-linear outer approximation** (`Za ≥ U·V`), where `U = Xmax−Xmin`, `V = Ymax−Ymin`.

---
## How to run

1. Install deps:

   ```bash
   pip install pyomo matplotlib
   ```

   You also need **GAMS with CPLEX** accessible on your PATH.

2. Run the script:

   ```bash
   python Main.py
   ```

---

## Model at a glance

### Sets

* `units`: list of plant units.

### Parameters

* `alpha[u], beta[u]` — sides (m) used in orientation swap.
* `Demin[i,j]` — **edge-to-edge** safety margin (m). (Symmetric.)
* Piping data per connected pair `(i,j)`:

  * `velocity[i,j]` (m/s), `Q[i,j]` (m³/s), `npp[i,j]` .
  * From these, a **per-meter cost** `c_per_m[i,j]` is computed.

### Decision variables

* `x[u], y[u]` — unit centres (m).
* `O[u] ∈ {0,1}` — orientation (swap `alpha`/`beta`).
* `l[u], d[u]` — side lengths implied by `O[u]`.
* `Xmin, Xmax, Ymin, Ymax` — land bounds (with `Xmin=Ymin=0` fixed).
* For each pair `(i,j)`:

  * `bx[i,j] ∈ {0,1}` — whether separation is enforced along X (else Y).
  * `sx[i,j], sy[i,j] ∈ {0,1}` — direction of separation along chosen axis.
  * `dx[i,j], dy[i,j] ≥ 0` — |xᵢ−xⱼ| and |yᵢ−yⱼ| (for piping length).
  * `t[i,j] = dx+dy` — L1 length used for pipe CAPEX.

### Constraints 

* **Orientation**: `l = alpha*O + beta*(1−O)`, `d = alpha+beta−l`.
* **Inside land**: each rectangle is within `[Xmin,Xmax]×[Ymin,Ymax]`.
* **Safety + non-overlap** (edge-to-edge): for every pair `(i,j)`
  If separated on X:
  `x[i] − x[j] ≥ 0.5(l[i]+l[j]) + Demin[i,j]` (or the symmetric form),
  else on Y with `d` replacing `l`.
  Big-M terms with binaries (`bx, sx, sy`) pick axis and direction.
* **Land area outer approx**: tangent cuts enforce `Za ≥ U·V`.

### Objective

```
min  LAND_PRICE * Za  +  Σ_{(i,j)∈pipes} c_per_m[i,j] * t[i,j]
```
**note** that piping is an important decision variable to consider in order to keep directly upstream buildings closer to eachother.
---

## Inputs you will edit most

### 1) Unit sizes

```python
alpha = {'compressor house':30, ...} #length of the building
beta  = {'compressor house':30, ...} #width of the building 
```

### 2) Safety distances (edge-to-edge)

Use `setD(i, j, value)` for each **unordered** pair you care about; the helper makes it symmetric:

```python
setD('compressor house','furnace',150) 
...
```
We calculated these as part of a research project in Safety and Loss Prevention with the use of PHAST (a specialised software)

### 3) Piping network

Use `set_pair` for `velocity` and `Q`, then set `npp`:

```python
set_pair(velocity, 'compressor house','furnace', 30.0) #3 or 30 depending on if the stream is liquid
set_pair(Q,        'compressor house','furnace', 116658.3/3600)
npp['compressor house','furnace'] = npp['furnace','compressor house'] = 1
```

**Note only pairs with both V, Q, and npp will translate to a full pipe in the problem**

---

## Cost models

* **Pipe per-meter cost** is derived from a diameter estimate (from `Q` & `velocity`) and escalated:

  ```
  PC_bare = BB * DIA^n2 * (CEPCI_2021/CEPCI_2006) * FX_rate
  PC_installed = PC_bare * (1 + F_install)
  c_per_m[i,j] = PC_installed * npp[i,j]
  ```
* **Land CAPEX** = `LAND_PRICE * area`, area approximated by `Za` which can be adapted base off of regional differences.

---

## Solver & limits

The script uses **Pyomo** with **GAMS/CPLEX**:

```python
solver = pyo.SolverFactory('gams:cplex')
solver.options['add_options'] = [
  'option reslim=300;',  # time limit (s)
  'option optcr=0.02;',  # relative optimality gap (2%)
]
res = solver.solve(m, tee=True)
```

Adjust `reslim` and `optcr` for faster or tighter runs.

---

## Customisation guide

* **Add a new pipe** between `A` and `B`:

  ```python
  set_pair(velocity, 'A','B', v)
  set_pair(Q,        'A','B', q)
  npp['A','B'] = npp['B','A'] = 1
  ```
* **Change a safety margin**:

  ```python
  setD('reactor','battery', 150)
  ```
* **Force a unit position** (optional anchors):

  ```python
  m.x['some unit'].fix(123.0)
  m.y['some unit'].fix(45.0)
  ```

---

## Citation 

Formulation adapted from common facility layout and process plant siting approaches using Pyomo
