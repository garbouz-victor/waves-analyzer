"""FEM edge-root contour reconstruction; no raster/marching image threshold."""
import numpy as np


def quadratic_edge_roots(a, mid, b, tol=1e-12):
    """Roots of the actual CG2 restriction to an edge, including all crossings."""
    # Endpoint values at a true mesh-vertex crossing can be ~3e-15 rather
    # than zero. Near tangency this produces two differently conditioned roots
    # on adjacent edges. Use the declared scalar root tolerance for endpoint
    # recognition, BEFORE the quadratic formula; no FEM coefficient is modified.
    if abs(a)<tol:
        a=0.0
    if abs(b)<tol:
        b=0.0
    c2,c1,c0=2*(a+b-2*mid),4*mid-3*a-b,a
    if max(abs(c2),abs(c1),abs(c0))<tol:
        return [0.0,1.0]  # identically zero edge
    roots=np.roots([c2,c1,c0] if abs(c2)>tol else [c1,c0])
    return sorted({float(np.clip(r.real,0,1)) for r in roots
                   if abs(r.imag)<tol and -tol<=r.real<=1+tol})


def contour_segments(solver, subdivisions=2):
    """Split reference triangles; roots are exact P1/P2, chords approximate arcs.

    Refine subdivisions independently to test geometric extraction. Return all
    segments; ambiguous cell topology is an error, never silently one crossing.
    """
    from dolfinx import fem
    import ufl
    subtriangles=[]
    n=subdivisions
    for i in range(n):
        for j in range(n-i):
            subtriangles.append(np.array([[i,j],[i+1,j],[i,j+1]],float)/n)
            if i+j<n-1:
                subtriangles.append(np.array([[i+1,j],[i+1,j+1],[i,j+1]],float)/n)
    refs=[]
    for tri in subtriangles:
        for k in range(3):
            a,b=tri[k],tri[(k+1)%3]
            refs.extend([a,(a+b)/2,b])
    refs=np.array(refs)
    mesh=solver.mesh
    cells=np.arange(mesh.topology.index_map(2).size_local,dtype=np.int32)
    vals=fem.Expression(ufl.split(solver.state)[2],refs).eval(mesh,cells).reshape(len(cells),len(refs))
    coords=mesh.geometry.x[mesh.geometry.dofmap[cells],:2]
    segments=[]
    for ci in range(len(cells)):
        # Only cells with sign change at subtriangle samples can contain a
        # resolved interface here. Tiny enclosed P2 loops require mesh refinement.
        for ti,tri in enumerate(subtriangles):
            values=vals[ci,ti*9:(ti+1)*9]
            if values.min()>0 or values.max()<0:
                continue
            crossings=[]
            zero_edges=[]
            for k in range(3):
                a,m,b=values[k*3:k*3+3]
                if max(abs(a),abs(m),abs(b))<1e-12:
                    zero_edges.append([tri[k],tri[(k+1)%3]])
                for root in quadratic_edge_roots(a,m,b):
                    r=tri[k]+root*(tri[(k+1)%3]-tri[k])
                    if not any(np.linalg.norm(r-v)<1e-10 for v in crossings):
                        crossings.append(r)
            candidates=zero_edges if zero_edges else ([crossings] if len(crossings)==2 else [])
            if not zero_edges and len(crossings)>2:
                raise RuntimeError(f"Ambiguous FEM contour in cell {ci}, subtriangle {ti}: "
                    f"roots={crossings}, edge samples={values.tolist()}; increase subdivisions or inspect tangency")
            for a,b in candidates:
                physical=[]
                for r in (a,b):
                    physical.append(coords[ci,0]+r[0]*(coords[ci,1]-coords[ci,0])+
                                    r[1]*(coords[ci,2]-coords[ci,0]))
                if np.linalg.norm(physical[0]-physical[1])>1e-12:
                    segments.append(physical)
    gathered=solver.comm.allgather(segments)
    return np.asarray([seg for rank in gathered for seg in rank],dtype=float).reshape(-1,2,2)


def connected_polylines(segments,tolerance=1e-9):
    """Return EVERY connected contour component, not just the first."""
    positions={}
    adjacency={}
    edges=set()
    for a,b in np.asarray(segments):
        keys=[tuple(np.rint(p/tolerance).astype(np.int64)) for p in (a,b)]
        if keys[0]==keys[1]:
            continue
        edge=tuple(sorted(keys))
        if edge in edges:
            continue
        edges.add(edge)
        for key,p in zip(keys,(a,b)):
            positions[key]=p
            adjacency.setdefault(key,set())
        adjacency[keys[0]].add(keys[1]);adjacency[keys[1]].add(keys[0])
    if any(len(v)>2 for v in adjacency.values()):
        raise RuntimeError("Branched/ambiguous contour connectivity")
    remaining=set(edges)
    lines=[]
    while remaining:
        involved={v for edge in remaining for v in edge}
        start=next((v for v in involved if len(adjacency[v])==1),next(iter(involved)))
        points=[positions[start]];current=start
        while True:
            candidates=[n for n in adjacency[current] if tuple(sorted((current,n))) in remaining]
            if not candidates:
                break
            nxt=candidates[0]
            remaining.remove(tuple(sorted((current,nxt))))
            points.append(positions[nxt]);current=nxt
        lines.append(np.array(points))
    return lines


def fit_circle(points):
    points=np.unique(np.asarray(points).reshape(-1,2),axis=0)
    if len(points)<6:
        raise ValueError("At least six distinct interface roots needed for circle fit")
    A=np.column_stack([2*points[:,0],2*points[:,1],np.ones(len(points))])
    xc,zc,c=np.linalg.lstsq(A,np.sum(points**2,axis=1),rcond=None)[0]
    radius=np.sqrt(c+xc*xc+zc*zc)
    residual=np.linalg.norm(points-[xc,zc],axis=1)-radius
    return {"center_x":float(xc),"center_z":float(zc),"radius":float(radius),
            "radial_rms":float(np.sqrt(np.mean(residual**2))),"points":len(points)}
