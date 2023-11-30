import sympy as sp

def christoffel_symbols(metric, coords):
    dims = metric.shape[0]
    christoffel = [[[0 for _ in range(dims)] for _ in range(dims)] for _ in range(dims)]
    inv_metric = metric.inv()
    for k in range(dims):
        for i in range(dims):
            for j in range(dims):
                sum_term = 0
                for l in range(dims):
                    partial_1 = sp.simplify(sp.diff(metric[j, l], coords[i]))
                    partial_2 = sp.simplify(sp.diff(metric[i, l], coords[j]))
                    partial_3 = sp.simplify(sp.diff(metric[i, j], coords[l]))
                    sum_term += inv_metric[k, l] * (partial_1 + partial_2 - partial_3)
                christoffel[k][i][j] = sp.simplify(sum_term / 2)
    return christoffel
def riemann_tensor(christoffel, coords):
    dims = len(coords)
    R = [[[sp.zeros(dims, dims) for _ in range(dims)] for _ in range(dims)] for _ in range(dims)]
    for rho in range(dims):
        for sigma in range(dims):
            for mu in range(dims):
                for nu in range(dims):
                    term1 = sp.simplify(sp.diff(christoffel[rho][nu][sigma], coords[mu]))
                    term2 = sp.simplify(sp.diff(christoffel[rho][mu][sigma], coords[nu]))
                    term3 = sp.simplify(sum(christoffel[rho][mu][lambda1] * christoffel[lambda1][nu][sigma] for lambda1 in range(dims)))
                    term4 = sp.simplify(sum(christoffel[rho][nu][lambda1] * christoffel[lambda1][mu][sigma] for lambda1 in range(dims)))
                    R[rho][sigma][mu][nu] = sp.simplify(term1 - term2 + term3 - term4)
    return R
def ricci_tensor(riemann, coords):
    dims = len(coords)
    ricci = sp.zeros(dims, dims)
    for mu in range(dims):
        for nu in range(dims):
            ricci[mu, nu] = sp.simplify(sum(riemann[lambda1][mu][lambda1][nu] for lambda1 in range(dims)))
    return ricci
def ricci_scalar(ricci, metric):
    inv_metric = metric.inv()
    dims = metric.shape[0]
    scalar = 0
    for mu in range(dims):
        for nu in range(dims):
            scalar += sp.simplify(inv_metric[mu, nu] * ricci[mu, nu])
    return sp.simplify(scalar)
def einstein_tensor(ricci_tensor, ricci_scalar, metric):
    dims = metric.shape[0]
    einstein = sp.zeros(dims, dims)
    for mu in range(dims):
        for nu in range(dims):
            einstein[mu, nu] = ricci_tensor[mu, nu] - sp.Rational(1/2) * metric[mu, nu] * ricci_scalar
    return sp.simplify(einstein)
def stress_energy_tensor_perfect_fluid(rho, p, u, metric):
    dims = metric.shape[0]
    stress_energy = sp.zeros(dims, dims)
    for mu in range(dims):
        for nu in range(dims):
            stress_energy[mu, nu] = (rho + p) * u[mu] * u[nu] + p * metric[mu, nu]
    return sp.simplify(stress_energy)
def einsteins_equations(einstein_tensor, stress_energy_tensor_perfect_fluid, metric):
    dims = metric.shape[0]
    einstein_eqs = sp.zeros(dims, dims)
    for mu in range(dims):
        for nu in range(dims):
            einstein_eqs[mu, nu] = einstein_tensor[mu, nu] - 8 * sp.pi * stress_energy_tensor_perfect_fluid[mu,nu]
    # Einstein's Equations: G_{\mu\nu} - 8πT_{\mu\nu} = 0 (In geometric units)
    # We return a tensor that should be all zeros if the equations are satisfied
    return sp.simplify(einstein_eqs)
