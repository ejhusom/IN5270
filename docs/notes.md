# IN5270 - Notes


## Finite element method


- W: Interval/domain. Divide this into N_e non-overlapping subintervals W(e),
  e=(0,...,N_e - 1)
- W(e): Each of these subdomains is called an *element*, identified by a unique
  number e.
- Nodes: On each element W(e) we introduce a set of points called *nodes*.
  Assume that the nodes are uniformly space throoughout the element, and that
  the boundary points of the elements are also nodes. Nodes are given numbers
  both within an element (local node numbers: r = (0, ..., d)) and in the
  global domain (global node numbers: i = (0, ..., N_n - 1)).
- Nodes and elements uniquely define a *finite element mesh*. A common special
  case is that of a uniformly partitioned mesh, where each element has the same
  length and the distance between nodes is constant.
- Finite element basis functions: phi_i(x), where the index corresponds to the
  global node number.
- We distinguish between internal nodes (in an element) and shared nodes.
