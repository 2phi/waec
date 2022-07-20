"""Fracture mechanics methods for WEak Layer AntiCrack nucleation model."""

# from scipy.optimize import root_scalar, basinhopping, brentq
# from scipy.integrate import quad, trapz, romberg


# def energy_criterion():
#     """Evaluate the energy criterion for a crack da = [xstart, xstop]."""
#     pass
#
#
# def stress_criterion(self, x, C, phi, li, ki, crit='quads', **kwargs):
#     """Evaluate the stress criterion at locations x."""
#     # Unused arguments
#     _ = kwargs
#
#     # Unpack strengths (given in kPa)
#     Sc = self.weak['Sc']
#     Tc = self.weak['Tc']
#
#     # Make sure x is np.ndarray and determine its size
#     x = np.asarray([x]) if np.isscalar(x) else np.asarray(x)
#     n = x.size
#
#     # Compute cumulate lengths list and segment index of all x
#     lic = np.cumsum(li)
#     nsegment = np.searchsorted(lic, x)
#
#     # Calculate global x coordinate of left segment ends and init stresses
#     x0 = np.insert(lic, 0, 0)
#     sig = np.zeros(n)
#     tau = np.zeros(n)
#
#     # Compute stresses at x and convert to kPa
#     for i, seg in enumerate(nsegment):
#         z0 = self.z(x[i] - x0[seg], C[:, [seg]], li[seg], phi, bed=ki[seg])
#         sig[i] = 1e3*max([0, self.sig(z0)])
#         tau[i] = 1e3*self.tau(z0)
#
#     # Evaluate stress criterion
#     if crit == 'quads':
#         return np.sqrt((sig/Sc)**2 + (tau/Tc)**2) - 1
#
# def external_potential(self, C, phi, **segments):
#     """Compute total external potential(PST or skier-on-slab setup)."""
#     # Rasterize solution
#     xq, zq, _ = self.rasterize_solution(C=C, phi=phi, **segments)
#     # Compute displacements where weight loads are applied
#     w0 = self.w(zq)
#     us = self.u(zq, z0=self.zs)
#     # Get weight loads
#     qn, qt = self.get_weight_load(phi)
#     # Integrate external work
#     Wext = trapz(qn*w0 + qt*us, xq)
#
#     if self.system == 'skier':
#         # Get skier weight and length of left free beam segment
#         m, l2 = segments['mi'][1], segments['li'][1]
#         # Compute solution at the skier's location
#         z0 = self.z(x=l2, C=C[:, [1]], l=l2, phi=phi, bed=False)
#         # Compute skier force
#         Fn, Ft = self.get_skier_load(m, phi)
#         # Add external work of the skier loading
#         Wext += Fn*self.w(z0) + Ft*self.u(z0, z0=-self.h/2)
#     elif self.system != 'pst-':
#         sys.exit('Input error: Only skier-on-slab and PST setups '
#                  + 'implemented at the moment.')
#
#     # Return potential of external forces Pext = - Wext
#     return -Wext
#
# def internal_potential(self, C, phi, **segments):
#     """Compute total internal potential(PST or skier-on-slab setup)."""
#     # Rasterize solution
#     xq, zq, xb = self.rasterize_solution(C=C, phi=phi, **segments)
#     # Compute section forces
#     N, M, V = self.N(zq), self.M(zq), self.V(zq)
#     # Compute stored energy of the slab (beam)
#     Pint = trapz(N**2/self.A11 + M**2/self.D11 + V**2/self.kA55, xq)/2
#
#     # Drop parts of the solution that are not a foundation
#     zweak = zq[:, ~np.isnan(xb)]
#     xweak = xb[~np.isnan(xb)]
#
#     if self.system == 'pst-':
#         # Compute displacments of segment on foundation
#         # w = self.w(zweak)
#         # u = self.u(zweak, z0=self.h/2)
#         eps = self.eps(zweak)
#         gamma = self.gamma(zweak)
#         sig = self.sig(zweak)
#         tau = self.tau(zweak)
#         # Compute stored energy of the weak layer (foundation)
#         # Pint += trapz(self.kn*w**2 + self.kt*u**2, xweak)/2
#         Pint += 1/2*trapz(sig*eps + tau*gamma, xweak)*self.t
#     elif self.system == 'skier':
#         # Split left and right bedded segments
#         zl, zr = np.array_split(zweak, 2, axis=1)
#         xl, xr = np.array_split(xweak, 2)
#         # Compute displacements on left and right foundations
#         wl, wr = self.w(zl), self.w(zr)
#         ul, ur = self.u(zl, z0=self.h/2), self.u(zr, z0=self.h/2)
#         # Compute stored energy of the weak layer (foundation)
#         Pint += trapz(self.kn*wl**2 + self.kt*ul**2, xl)/2
#         Pint += trapz(self.kn*wr**2 + self.kt*ur**2, xr)/2
#     else:
#         sys.exit('Input error: Only skier-on-slab and PST setups '
#                  + 'implemented at the moment.')
#
#     return Pint

# def find_roots(self, L, C, phi, nintervals=50, **kwargs):
#     """Find points where the stress criterion is satisfied identically."""
#     # Subdivide domain into intervals
#     intvls = np.linspace(0, L, num=nintervals+1)
#     roots = []
#
#     # See if we can find roots in given intervals
#     for i in range(nintervals):
#         try:
#             roots.append(root_scalar(
#                 self.stress_criterion, method='brentq',
#                 args=(C, phi, kwargs['li'], kwargs['ki']),
#                 bracket=intvls[i:i+2], xtol=1e-1).root)
#         except ValueError:
#             pass
#
#     # Determine index to insert skier position
#     iskier = np.searchsorted(roots, L/2)
#     # Add domain ends and skier position to points of interest
#     poi = np.unique(np.insert(roots, [0, iskier, len(roots)], [0, L/2, L]))
#     # Compute new segment lengths
#     li = poi[1:] - poi[:-1]
#     # Compute new segment midpoints
#     xmid = (poi[1:] + poi[:-1])/2
#     # Compute new foundation order
#     ki = self.stress_criterion(
#         xmid, C, phi, kwargs['li'], kwargs['ki']) < 0
#     k0 = np.full_like(ki, True)
#     # Compute new position of skier weight
#     lic = np.cumsum(li)[:-1]
#     iskier = np.searchsorted(lic, L/2)
#     mi = np.insert(np.zeros_like(lic), iskier, kwargs['mi'].sum())
#     print(lic)
#     print(mi)
#
#     # Assemble update segmentation as output
#     segments = {
#         'nocrack': {'li': li, 'mi': mi, 'ki': k0},
#         'crack': {'li': li, 'mi': mi, 'ki': ki},
#         'both': {'li': li, 'mi': mi, 'ki': ki, 'k0': k0}}
#
#     return segments
