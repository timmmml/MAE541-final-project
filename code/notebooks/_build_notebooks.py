"""Build the 00–10 notebook series from compact definitions.

Run this script (from anywhere) to (re)generate all ``code/notebooks/*.ipynb``
files. Each notebook is defined below as a list of ``("md", text)`` /
``("code", text)`` tuples; ``nbformat`` assembles them into an ``.ipynb`` JSON
file. This keeps the source-of-truth concise and avoids hand-editing JSON.
"""

from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import nbformat as nbf


# ---------------------------------------------------------------------------
# Common preamble used by every notebook.
# ---------------------------------------------------------------------------

PREAMBLE = dedent(
    """\
    %matplotlib inline
    import numpy as np
    import matplotlib.pyplot as plt

    from steering.params import ModelParams, ForcingParams
    from steering.models import (
        DuffingModel,
        BesselSteeringModel,
        ContinuousPFLModel,
        DiscretePFLModel,
        FullCircuitModel,
    )
    from steering.dynamics import VelocityDynamics, AccelerationDynamics
    from steering.integrator import Simulation
    from steering.visualization.style import use_paper_style

    use_paper_style()
    """
)


# ---------------------------------------------------------------------------
# Notebook content definitions. Each is a list of (kind, src) tuples.
# kind = 'md' or 'code'.
# ---------------------------------------------------------------------------


def nb_00():
    return [
        ("md", "# 00 — Model Hierarchy Tour\n\n"
               "Verify the codebase loads, instantiate every model at a common "
               "parameter set, and confirm the §9.3 consistency identities "
               "hold numerically."),
        ("code", PREAMBLE),
        ("md", "## Common parameters\n\n"
               "We pick $\\kappa_h = \\kappa_g = 3$, $\\delta = \\pi/4$, "
               "$S = A = 1$, $\\Delta = 3\\pi/8$ (the *Drosophila* PFL3 "
               "phase shift)."),
        ("code", dedent("""\
            params = ModelParams(kappa_h=3.0, kappa_g=3.0, delta=np.pi/4)
            params_disc8 = params.replace(N_neurons=8)
            params_disc64 = params.replace(N_neurons=64)

            duffing = DuffingModel.from_params(params)
            bessel = BesselSteeringModel()
            cont = ContinuousPFLModel(n_quad=128)
            disc = DiscretePFLModel()
            full = FullCircuitModel()
            print('c1, c3 =', bessel.taylor_coefficients(params))
        """)),
        ("md", "## $U(\\theta)$ across all five models\n\n"
               "Bessel and ContinuousPFL with the quadratic nonlinearity "
               "should overlap exactly. Duffing should match near $\\theta=0$ "
               "and deviate at the boundaries. DiscretePFL with $N=64$ should "
               "be visually indistinguishable from the continuous model; "
               "$N=8$ shows discretization artifacts."),
        ("code", dedent("""\
            theta = np.linspace(-np.pi, np.pi, 401)
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(theta, bessel.steering_drive(theta, params), label='Bessel', lw=2)
            ax.plot(theta, cont.steering_drive(theta, params), label='ContinuousPFL', ls=':', lw=2)
            ax.plot(theta, disc.steering_drive(theta, params_disc64), label='DiscretePFL N=64', ls='--')
            ax.plot(theta, disc.steering_drive(theta, params_disc8), label='DiscretePFL N=8', alpha=0.7)
            ax.plot(theta, duffing.steering_drive(theta), label='Duffing (cubic)', ls='-.', lw=1.2)
            ax.axhline(0, color='0.6', lw=0.5)
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$U(\\theta)$')
            ax.legend(); plt.show()
        """)),
        ("md", "## Consistency table (§9.3 identity chains)"),
        ("code", dedent("""\
            grid = np.linspace(-np.pi+0.05, np.pi-0.05, 21)
            Ub = np.asarray(bessel.steering_drive(grid, params))
            Uc = np.asarray(cont.steering_drive(grid, params))
            Ud = duffing.steering_drive(grid)
            print(f'max|Bessel - Continuous|        = {np.max(np.abs(Ub-Uc)):.2e}')
            print(f'max|Duffing - Bessel| / |Bessel| = {np.max(np.abs(Ud-Ub))/np.max(np.abs(Ub)):.2e}')
            Ud256 = disc.steering_drive(grid, params.replace(N_neurons=256))
            print(f'max|Continuous - Discrete N=256| = {np.max(np.abs(Uc-Ud256)):.2e}')
            full_w0 = full.steering_drive(grid, params)
            print(f'max|FullCircuit(W=0) - Continuous| = {np.max(np.abs(full_w0-Uc)):.2e}')
        """)),
        ("md", "## Effect of nonlinearity\n\n"
               "Pick a bistable parameter set; compare $U(\\theta)$ for "
               "$f \\in \\{x^2, \\mathrm{ELU}, \\mathrm{ReLU}, "
               "\\mathrm{softplus}\\}$ using the continuous model. The "
               "Bessel closed form is *only* exact for $f=x^2$."),
        ("code", dedent("""\
            p_bi = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            theta = np.linspace(-np.pi, np.pi, 401)
            fig, ax = plt.subplots(figsize=(7, 4))
            for nl in ['quadratic', 'elu', 'relu', 'softplus']:
                U = ContinuousPFLModel(n_quad=128).steering_drive(theta, p_bi.replace(nonlinearity=nl))
                ax.plot(theta, U / np.max(np.abs(U) + 1e-12), label=nl)
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$U(\\theta)$ (normalised)')
            ax.legend(); ax.axhline(0, color='0.6', lw=0.5)
            plt.show()
        """)),
    ]


def nb_01():
    return [
        ("md", "# 01 — Steering Profiles $U(\\theta)$ and $V(\\theta)$\n\n"
               "Build geometric intuition: scan $(\\theta, \\delta)$ "
               "heatmaps, watch a single well split into a double well as "
               "$\\delta$ crosses the pitchfork, and overlay Duffing on "
               "Bessel to mark where the cubic approximation breaks down."),
        ("code", PREAMBLE),
        ("md", "## $U(\\theta, \\delta)$ heatmap (Bessel)\n\n"
               "Fix $\\kappa = 2$ and scan $(\\theta, \\delta)$ on a grid. "
               "The zero contours of this heatmap *are* the fixed-point "
               "branches as $\\delta$ varies — a direct visualisation of "
               "the pitchfork."),
        ("code", dedent("""\
            kappa = 2.0
            thetas = np.linspace(-np.pi, np.pi, 251)
            deltas = np.linspace(0.0, np.pi/2, 121)
            U = np.empty((deltas.size, thetas.size))
            bessel = BesselSteeringModel()
            for i, d in enumerate(deltas):
                p = ModelParams(kappa_h=kappa, kappa_g=kappa, delta=d)
                U[i] = bessel.steering_drive(thetas, p)
            fig, ax = plt.subplots(figsize=(7, 4.5))
            mesh = ax.pcolormesh(thetas, deltas, U, cmap='RdBu_r', shading='auto', vmin=-np.max(np.abs(U)), vmax=np.max(np.abs(U)))
            ax.contour(thetas, deltas, U, levels=[0.0], colors='k', linewidths=1.0)
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$\\delta$')
            fig.colorbar(mesh, ax=ax, label=r'$U(\\theta,\\delta)$')
            plt.show()
        """)),
        ("md", "## $V(\\theta)$ cross-sections through the pitchfork\n\n"
               "Pick a few $\\delta$ values: below, at, above the bifurcation. "
               "Watch the single well split into two."),
        ("code", dedent("""\
            chosen = [0.5, 1.20, 1.34, 1.45]
            thetas = np.linspace(-2.0, 2.0, 401)
            fig, ax = plt.subplots(figsize=(7, 4))
            for d in chosen:
                p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=d)
                V = bessel.steering_potential(thetas, p)
                ax.plot(thetas, V, label=fr'$\\delta={d}$')
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$V(\\theta)$')
            ax.legend(); plt.show()
        """)),
        ("md", "## Bessel vs Duffing overlay\n\n"
               "On each cross-section, overlay the Duffing cubic. The error "
               "grows with $|\\theta|$ — quantify the relative L^infinity "
               "error in a band around the origin."),
        ("code", dedent("""\
            fig, ax = plt.subplots(figsize=(7, 4))
            for d in [1.34, 1.40, 1.50]:
                p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=d)
                duff = DuffingModel.from_params(p)
                Ub = bessel.steering_drive(thetas, p)
                Ud = duff.steering_drive(thetas)
                ax.plot(thetas, Ub, label=fr'Bessel $\\delta={d}$')
                ax.plot(thetas, Ud, ls='--', alpha=0.7)
            ax.axhline(0, color='0.6', lw=0.5)
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$U(\\theta)$')
            ax.legend(); plt.show()
        """)),
        ("md", "## Sharpening bumps: $\\kappa$ at fixed $\\delta$\n\n"
               "Higher $\\kappa$ deepens and narrows the wells. Verify by "
               "plotting $V(\\theta)$ for several $\\kappa$ values."),
        ("code", dedent("""\
            fig, ax = plt.subplots(figsize=(7, 4))
            for kappa in [1.0, 2.0, 3.0, 5.0]:
                p = ModelParams(kappa_h=kappa, kappa_g=kappa, delta=1.45)
                V = bessel.steering_potential(thetas, p)
                ax.plot(thetas, V / np.max(np.abs(V) + 1e-12), label=fr'$\\kappa={kappa}$')
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$V/|V|_\\infty$')
            ax.legend(); plt.show()
        """)),
    ]


def nb_02():
    return [
        ("md", "# 02 — Pitchfork Bifurcation Analysis\n\n"
               "Compute the analytical $c_1, c_3$ Taylor coefficients on a "
               "$(\\kappa, \\delta)$ grid; locate the pitchfork as the "
               "$c_1=0$ contour; classify supercritical vs subcritical via "
               "$\\mathrm{sgn}(c_3)$; overlay the numerical fixed-point "
               "branches; verify the broad-bump anchor "
               "$\\sin\\Delta\\cos\\delta = 0$ in the $\\kappa\\to 0$ limit."),
        ("code", PREAMBLE),
        ("code", dedent("""\
            from steering.analysis.bifurcation import (
                pitchfork_bifurcation_curve,
                numerical_bifurcation_diagram_1d,
            )
            from steering.visualization.bifurcation_plot import plot_bifurcation_2d

            bessel = BesselSteeringModel()
        """)),
        ("md", "## Bifurcation diagram on $(\\kappa, \\delta)$\n\n"
               "We sweep $\\kappa_h$ over $(\\kappa_g$ slaved to it via "
               "post-processing), and $\\delta$ over $[0, \\pi/2]$."),
        ("code", dedent("""\
            kappas = np.linspace(0.5, 6.0, 61)
            deltas = np.linspace(0.0, np.pi/2 - 0.001, 61)
            base = ModelParams(kappa_h=1.0, kappa_g=1.0, delta=0.0)
            c1_grid = np.empty((kappas.size, deltas.size))
            c3_grid = np.empty_like(c1_grid)
            for i, k in enumerate(kappas):
                for j, d in enumerate(deltas):
                    p = base.replace(kappa_h=k, kappa_g=k, delta=d)
                    c1, c3 = bessel.taylor_coefficients(p)
                    c1_grid[i, j] = c1
                    c3_grid[i, j] = c3

            from steering.analysis.bifurcation import BifurcationDiagram
            diag = BifurcationDiagram(
                'kappa', 'delta', kappas, deltas, c1_grid, c3_grid,
            )
            fig, ax = plt.subplots(figsize=(7, 4.5))
            plot_bifurcation_2d(diag, r'$\\kappa$', r'$\\delta$', ax=ax)
            plt.show()
        """)),
        ("md", "## Broad-bump limit anchor\n\n"
               "For $\\kappa \\to 0$, the bifurcation condition reduces to "
               "$\\sin\\Delta\\cos\\delta = 0$, i.e. $\\delta = \\pi/2$ (with "
               "$\\Delta = 3\\pi/8 \\neq 0$). Compute and verify."),
        ("code", dedent("""\
            for k in [0.05, 0.1, 0.5, 1.0, 2.0, 5.0]:
                # Find delta where c1 = 0 by bisection.
                from scipy.optimize import brentq
                def c1_of(d, k=k):
                    return bessel.taylor_coefficients(base.replace(kappa_h=k, kappa_g=k, delta=d))[0]
                try:
                    d_star = brentq(c1_of, 0.5, np.pi/2 - 1e-3)
                except ValueError:
                    d_star = float('nan')
                print(f'kappa = {k:.2f}: pitchfork at delta = {d_star:.4f} (pi/2 = {np.pi/2:.4f})')
        """)),
        ("md", "## 1D bifurcation: stable/unstable branches as $\\delta$ varies\n\n"
               "Fix $\\kappa = 2$ and find all fixed points of $U(\\theta) = 0$ "
               "for each $\\delta$. Show stable branches as solid points, "
               "unstable as open."),
        ("code", dedent("""\
            data = numerical_bifurcation_diagram_1d(
                bessel, 'delta', np.linspace(0.5, np.pi/2 - 0.01, 80),
                base.replace(kappa_h=2.0, kappa_g=2.0),
                theta_range=(-np.pi, np.pi),
            )
            fig, ax = plt.subplots(figsize=(7, 4))
            for d, ths, stabs in zip(data['param'], data['thetas'], data['stabilities']):
                for th, s in zip(ths, stabs):
                    color = 'C0' if s == 'stable' else 'C3'
                    ax.plot(d, th, '.', color=color, ms=3)
            ax.set_xlabel(r'$\\delta$'); ax.set_ylabel(r'$\\theta^*$')
            ax.set_title('1D bifurcation: stable (blue) vs unstable (red)')
            plt.show()
        """)),
        ("md", "## Critical damping for the acceleration system\n\n"
               "$\\gamma_c = 2\\sqrt{|c_1|}$ at any fixed point. Above this, "
               "approach is overdamped; below, oscillatory."),
        ("code", dedent("""\
            for d in [1.34, 1.40, 1.45, 1.50]:
                p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=d)
                c1, c3 = bessel.taylor_coefficients(p)
                if c1 > 0 and c3 < 0:
                    # Stability eigenvalues at well: U' = c1 + 3 c3 theta_*^2 = -2 c1 (the well at sqrt(-c1/c3)).
                    Up_at_well = -2 * c1
                    gamma_c = 2 * np.sqrt(abs(Up_at_well))
                    print(f'delta = {d:.2f}: gamma_c = {gamma_c:.3f}')
        """)),
    ]


def nb_03():
    return [
        ("md", "# 03 — Phase Portraits and Hamiltonian Structure\n\n"
               "Visualise the 2D phase space of the acceleration system. "
               "Compare Duffing vs Bessel in monostable, weakly bistable, "
               "and strongly bistable regimes."),
        ("code", PREAMBLE),
        ("code", dedent("""\
            from steering.visualization.phase_portrait import plot_phase_portrait
            from steering.analysis.fixed_points import find_fixed_points_2d

            params_sets = [
                ('monostable',  ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.0)),
                ('weak bistable', ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.36)),
                ('strong bistable', ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.55)),
            ]
        """)),
        ("md", "## Print c1, c3, theta_max for each regime"),
        ("code", dedent("""\
            bessel = BesselSteeringModel()
            for name, p in params_sets:
                c1, c3 = bessel.taylor_coefficients(p)
                tmax = np.sqrt(2 * c1 / abs(c3)) if (c1 > 0 and c3 < 0) else float('nan')
                print(f'{name:<18}: c1={c1:+.3f}  c3={c3:+.3f}  theta_max={tmax:.3f}')
        """)),
        ("md", "## Phase portraits: Duffing vs Bessel, undamped\n\n"
               "Hamiltonian level sets in light blue; saddle separatrices in "
               "black; fixed points marked."),
        ("code", dedent("""\
            fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
            for j, (name, p) in enumerate(params_sets):
                duff = DuffingModel.from_params(p)
                for i, (label, model) in enumerate([('Duffing', duff), ('Bessel', bessel)]):
                    dyn = AccelerationDynamics(model=model, gamma=0.0, topology='planar')
                    plot_phase_portrait(
                        dyn, p, theta_range=(-np.pi, np.pi), v_range=(-6, 6),
                        n_grid=18, n_contours=12, ax=axes[i, j],
                    )
                    axes[i, j].set_title(f'{label}, {name}')
            fig.tight_layout(); plt.show()
        """)),
        ("md", "## Damped trajectories\n\n"
               "Add $\\gamma = 0.1$ and overlay several initial-condition "
               "trajectories — they spiral into the stable foci."),
        ("code", dedent("""\
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            dyn = AccelerationDynamics(model=bessel, gamma=0.1)
            sim = Simulation(dyn, p)
            fig, ax = plt.subplots(figsize=(6, 5))
            plot_phase_portrait(dyn, p, theta_range=(-np.pi, np.pi),
                                v_range=(-5, 5), n_grid=18, n_contours=10,
                                overlay_separatrix=False, ax=ax)
            for ic in [(-1.5, 1.0), (1.5, -1.0), (-2.5, 0.0), (0.5, 2.5)]:
                res = sim.run(np.array(ic), (0.0, 60.0))
                ax.plot(res.states[:, 0], res.states[:, 1], 'C3-', lw=0.9, alpha=0.6)
            ax.set_xlim(-np.pi, np.pi); ax.set_ylim(-5, 5)
            plt.show()
        """)),
    ]


def nb_04():
    return [
        ("md", "# 04 — Homoclinic Orbits\n\n"
               "Construct the analytical Duffing sech homoclinic and the "
               "numerical Bessel homoclinic; compare them; quantify the "
               "Duffing-validity boundary."),
        ("code", PREAMBLE),
        ("code", dedent("""\
            from steering.analysis.homoclinic import numerical_homoclinic
            bessel = BesselSteeringModel()
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            c1, c3 = bessel.taylor_coefficients(p)
            duff = DuffingModel(c1=c1, c3=c3)
            print(f'c1={c1:.3f}, c3={c3:.3f}, theta_max(Duffing)={(2*c1/abs(c3))**0.5:.3f}')
        """)),
        ("md", "## Analytical Duffing homoclinic"),
        ("code", dedent("""\
            t = np.linspace(-6, 6, 1201)
            th_d, v_d = duff.homoclinic_orbit(t)
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(t, th_d, label=r'$\\theta_0(t)$')
            ax.plot(t, v_d, label=r'$v_0(t)$')
            ax.legend(); ax.set_xlabel('t'); plt.show()
        """)),
        ("md", "## Numerical Bessel homoclinic (apex-mirror method)"),
        ("code", dedent("""\
            t_b, th_b, v_b = numerical_homoclinic(bessel, p)
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(t_b, th_b, label='Bessel theta')
            th_d2, _ = duff.homoclinic_orbit(t_b)
            ax.plot(t_b, th_d2, label='Duffing theta', ls='--')
            ax.legend(); ax.set_xlabel('t'); plt.show()
        """)),
        ("md", "## Pointwise error along the orbit"),
        ("code", dedent("""\
            U_b = np.asarray(bessel.steering_drive(th_b, p))
            U_d = duff.steering_drive(th_b)
            fig, ax = plt.subplots(figsize=(6, 3))
            ax.semilogy(t_b, np.abs(U_b - U_d) + 1e-15)
            ax.set_xlabel('t')
            ax.set_ylabel(r'$|U_{\\mathrm{Bessel}} - U_{\\mathrm{Duffing}}|$ on the orbit')
            plt.show()
        """)),
        ("md", "## Validity boundary: relative L^2 error vs $\\theta_{\\max}/\\pi$"),
        ("code", dedent("""\
            from scipy.integrate import simpson
            deltas = np.linspace(1.34, 1.55, 12)
            errs, tmaxes = [], []
            for d in deltas:
                pp = p.replace(delta=d)
                c1d, c3d = bessel.taylor_coefficients(pp)
                if not (c1d > 0 and c3d < 0):
                    continue
                duf = DuffingModel(c1=c1d, c3=c3d)
                t_n, th_n, _ = numerical_homoclinic(bessel, pp)
                th_anal, _ = duf.homoclinic_orbit(t_n)
                num = simpson((th_n - th_anal)**2, x=t_n)
                den = simpson(th_n**2, x=t_n)
                errs.append(np.sqrt(num/den))
                tmaxes.append(np.max(np.abs(th_n)) / np.pi)
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(tmaxes, errs, 'o-')
            ax.set_xlabel(r'$\\theta_{\\max}/\\pi$')
            ax.set_ylabel(r'rel. L^2 error (Duffing vs Bessel orbit)')
            plt.show()
        """)),
    ]


# ---------------------------------------------------------------------------
# Notebooks 05-10 below.
# ---------------------------------------------------------------------------


def nb_05():
    return [
        ("md", "# 05 — Melnikov Function and Chaos Threshold\n\n"
               "Plot $M(t_0)$ analytically and numerically; compute "
               "$F_{\\mathrm{crit}}$ vs $\\omega$ and vs $\\delta$; validate "
               "by direct simulation. Note the spec's $F_{\\mathrm{crit}}$ "
               "formula is corrected here — the leading factor is "
               "$c_1^{3/2}$ (not $c_1$)."),
        ("code", PREAMBLE),
        ("code", dedent("""\
            from steering.analysis.melnikov import (
                melnikov_analytical,
                melnikov_critical_forcing,
                melnikov_critical_forcing_numerical,
                melnikov_numerical,
            )
            from steering.analysis.homoclinic import numerical_homoclinic

            bessel = BesselSteeringModel()
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            c1, c3 = bessel.taylor_coefficients(p)
            print(f'c1={c1:.3f}, c3={c3:.3f}')
        """)),
        ("md", "## Analytical $M(t_0)$ for several $(F, \\omega)$"),
        ("code", dedent("""\
            gamma = 0.05
            t0 = np.linspace(0, 4*np.pi, 401)
            fig, ax = plt.subplots(figsize=(7, 4))
            for F in [0.05, 0.10, 0.20]:
                for omega in [1.0]:
                    M = melnikov_analytical(c1, c3, gamma, omega, F, t0)
                    ax.plot(t0, M, label=fr'$F={F}, \\omega={omega}$')
            ax.axhline(0, color='0.5', lw=0.5)
            ax.set_xlabel(r'$t_0$'); ax.set_ylabel(r'$M(t_0)$')
            ax.legend(); plt.show()
        """)),
        ("md", "## $F_{\\mathrm{crit}}(\\omega)$ — Duffing vs Bessel"),
        ("code", dedent("""\
            omegas = np.linspace(0.2, 4.0, 30)
            F_duff = np.array([melnikov_critical_forcing(c1, c3, gamma, om) for om in omegas])
            homo = numerical_homoclinic(bessel, p)
            F_bess = np.array([melnikov_critical_forcing_numerical(bessel, p, gamma, om, homoclinic=homo) for om in omegas])
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.semilogy(omegas, F_duff, label='Duffing analytical')
            ax.semilogy(omegas, F_bess, '--', label='Bessel numerical')
            ax.set_xlabel(r'$\\omega$'); ax.set_ylabel(r'$F_{\\mathrm{crit}}$')
            ax.legend(); plt.show()
        """)),
        ("md", "## $F_{\\mathrm{crit}}(\\delta)$ at fixed $\\omega$\n\n"
               "Across $\\delta$ values where the system is bistable. The "
               "minimum is the most chaos-susceptible separation."),
        ("code", dedent("""\
            deltas = np.linspace(1.33, 1.55, 18)
            F_duff = []
            for d in deltas:
                pp = p.replace(delta=d)
                c1d, c3d = bessel.taylor_coefficients(pp)
                if c1d > 0 and c3d < 0:
                    F_duff.append(melnikov_critical_forcing(c1d, c3d, gamma, 1.0))
                else:
                    F_duff.append(np.nan)
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(deltas, F_duff)
            ax.set_xlabel(r'$\\delta$'); ax.set_ylabel(r'$F_{\\mathrm{crit}}$')
            plt.show()
        """)),
        ("md", "## Direct simulation validation\n\n"
               "Pick $\\omega=1$, $\\gamma=0.05$, and three $F$ levels: "
               "below, near, and above $F_{\\mathrm{crit}}$. Run long sims "
               "and inspect the Poincaré section."),
        ("code", dedent("""\
            from steering.analysis.poincare import stroboscopic_section
            F_crit = melnikov_critical_forcing(c1, c3, gamma, 1.0)
            print(f'F_crit (Duffing) = {F_crit:.4f}')
            dyn = AccelerationDynamics(model=bessel, gamma=gamma)
            fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
            for ax, F in zip(axes, [0.5*F_crit, 1.5*F_crit, 4.0*F_crit]):
                forc = ForcingParams(F=F, omega=1.0)
                sim = Simulation(dyn, p, forc, rtol=1e-9, atol=1e-11)
                T = 2*np.pi
                res = sim.run(np.array([0.5, 0.0]), (0.0, 800*T), dense_output=True)
                sec = stroboscopic_section(res, omega=1.0, transient_periods=200)
                ax.plot(sec[:, 0], sec[:, 1], 'k.', ms=1.0, alpha=0.5)
                ax.set_title(fr'$F={F:.3f}$')
                ax.set_xlabel(r'$\\theta$')
            axes[0].set_ylabel(r'$v$'); plt.show()
        """)),
    ]


def nb_06():
    return [
        ("md", "# 06 — Poincaré Sections and Route to Chaos\n\n"
               "Stroboscopic Poincaré bifurcation diagrams; visualize the "
               "strange attractor; compute the maximal Lyapunov exponent."),
        ("code", PREAMBLE),
        ("code", dedent("""\
            from steering.analysis.poincare import bifurcation_diagram_poincare
            from steering.analysis.lyapunov import max_lyapunov_exponent

            bessel = BesselSteeringModel()
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            gamma, omega = 0.1, 1.0
        """)),
        ("md", "## Bifurcation diagram: stroboscopic $\\theta$ vs $F$"),
        ("code", dedent("""\
            dyn = AccelerationDynamics(model=bessel, gamma=gamma)
            F_vals = np.linspace(0.05, 1.6, 50)
            diag = bifurcation_diagram_poincare(
                dyn, bessel, p, gamma=gamma, omega=omega,
                sweep_param='F', sweep_values=F_vals,
                state0=np.array([0.5, 0.0]),
                n_periods_transient=120, n_periods_record=60,
            )
            fig, ax = plt.subplots(figsize=(8, 4.5))
            for v, vals in zip(diag['sweep_values'], diag['thetas']):
                ax.plot(np.full_like(vals, v), vals, 'k.', ms=0.7, alpha=0.5)
            ax.set_xlabel(r'$F$'); ax.set_ylabel(r'$\\theta_{\\mathrm{strobe}}$')
            plt.show()
        """)),
        ("md", "## Lyapunov exponent across the same sweep"),
        ("code", dedent("""\
            F_lyap = np.linspace(0.05, 1.6, 12)
            lambdas = []
            for F in F_lyap:
                lam = max_lyapunov_exponent(
                    dyn, p, gamma=gamma, F=F, omega=omega,
                    state0=np.array([0.5, 0.0]),
                    t_total=200.0, renorm_interval=1.0,
                )
                lambdas.append(lam)
            fig, ax = plt.subplots(figsize=(7, 3.5))
            ax.plot(F_lyap, lambdas, 'o-')
            ax.axhline(0, color='0.5', lw=0.5)
            ax.set_xlabel(r'$F$'); ax.set_ylabel(r'$\\lambda_{\\max}$')
            plt.show()
        """)),
    ]


def nb_07():
    return [
        ("md", "# 07 — Topology Comparison: $\\mathbb{R}^2$ vs $S^1 \\times \\mathbb{R}$\n\n"
               "When does the planar approximation break down and the "
               "cylindrical topology produce qualitatively different "
               "dynamics?"),
        ("code", PREAMBLE),
        ("code", dedent("""\
            bessel = BesselSteeringModel()
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
        """)),
        ("md", "## Trajectory comparison\n\n"
               "Strongly bistable + strong forcing → trajectory may wander "
               "to neighbouring wells via the $\\theta = \\pi$ saddle on "
               "$S^1$. Inside the Duffing-valid region the two topologies "
               "agree."),
        ("code", dedent("""\
            forc = ForcingParams(F=0.6, omega=1.0)
            dyn_planar = AccelerationDynamics(model=bessel, gamma=0.05, topology='planar')
            dyn_cylin = AccelerationDynamics(model=bessel, gamma=0.05, topology='cylindrical')
            sim_p = Simulation(dyn_planar, p, forc, rtol=1e-9, atol=1e-11)
            sim_c = Simulation(dyn_cylin, p, forc, rtol=1e-9, atol=1e-11)
            r_p = sim_p.run(np.array([0.5, 0.0]), (0.0, 200.0))
            r_c = sim_c.run(np.array([0.5, 0.0]), (0.0, 200.0))
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.plot(r_p.t, r_p.states[:, 0], label='planar')
            ax.plot(r_c.t, r_c.states[:, 0], label='cylindrical', alpha=0.7)
            ax.set_xlabel('t'); ax.set_ylabel(r'$\\theta$'); ax.legend()
            plt.show()
        """)),
        ("md", "## Winding number on $S^1$"),
        ("code", dedent("""\
            # Run on cylindrical topology but accumulate the *unwrapped* theta
            # by integrating planar; then compute (unwrapped(t) - unwrapped(0))/(2pi).
            r_long = sim_p.run(np.array([0.5, 0.0]), (0.0, 1000.0))
            theta_unw = r_long.states[:, 0]
            W_t = (theta_unw - theta_unw[0]) / (2.0*np.pi)
            fig, ax = plt.subplots(figsize=(7, 3))
            ax.plot(r_long.t, W_t)
            ax.set_xlabel('t'); ax.set_ylabel('winding number $W(t)$')
            plt.show()
        """)),
    ]


def nb_08():
    return [
        ("md", "# 08 — Finite-$N$ (Discrete) Effects\n\n"
               "Convergence of the discrete population sum to the continuous "
               "integral; bifurcation locus shift at biological $N=8$; "
               "symmetry preservation at finite $N$."),
        ("code", PREAMBLE),
        ("code", dedent("""\
            cont = ContinuousPFLModel(n_quad=256)
            disc = DiscretePFLModel()
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            theta = np.linspace(-np.pi, np.pi, 401)
            U_inf = cont.steering_drive(theta, p)
        """)),
        ("md", "## $L^\\infty$ convergence of $U(\\theta)$"),
        ("code", dedent("""\
            Ns = [4, 8, 16, 32, 64, 256]
            errs = []
            for N in Ns:
                U_N = disc.steering_drive(theta, p.replace(N_neurons=N))
                errs.append(np.max(np.abs(U_N - U_inf)))
            fig, ax = plt.subplots(figsize=(6, 3.5))
            ax.loglog(Ns, errs, 'o-')
            ax.set_xlabel('N'); ax.set_ylabel(r'$\\|U_N - U_\\infty\\|_\\infty$')
            plt.show()
        """)),
        ("md", "## Bifurcation locus at finite $N$"),
        ("code", dedent("""\
            from scipy.optimize import brentq
            from steering.models import BesselSteeringModel
            bessel = BesselSteeringModel()

            def find_pitchfork(model, kappa, get_U_at_theta):
                def c1_at(d):
                    pp = ModelParams(kappa_h=kappa, kappa_g=kappa, delta=d)
                    h = 1e-3
                    return (get_U_at_theta(h, pp) - get_U_at_theta(-h, pp)) / (2*h)
                try:
                    return brentq(c1_at, 0.5, np.pi/2 - 1e-3)
                except ValueError:
                    return float('nan')

            kappas = np.linspace(0.5, 5, 12)
            curves = {}
            for N in [None, 8, 16, 64]:
                ds = []
                for k in kappas:
                    if N is None:
                        ds.append(find_pitchfork(bessel, k,
                            lambda th, pp: float(np.asarray(bessel.steering_drive(th, pp)))))
                    else:
                        ds.append(find_pitchfork(disc, k,
                            lambda th, pp, _N=N: float(disc.steering_drive(th, pp.replace(N_neurons=_N)))))
                curves[f'N={N}' if N is not None else 'continuous'] = ds
            fig, ax = plt.subplots(figsize=(6, 4))
            for label, ds in curves.items():
                ax.plot(kappas, ds, label=label)
            ax.set_xlabel(r'$\\kappa$'); ax.set_ylabel(r'$\\delta_{\\mathrm{pitchfork}}$')
            ax.legend(); plt.show()
        """)),
    ]


def nb_09():
    return [
        ("md", "# 09 — Full Circuit with PFL2 Modulation\n\n"
               "How does the indirect pathway (PFL2 → DNa03) modulate the "
               "steering drive, the bifurcation, and the chaos threshold?"),
        ("code", PREAMBLE),
        ("code", dedent("""\
            full = FullCircuitModel()
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            theta = np.linspace(-np.pi, np.pi, 401)
        """)),
        ("md", "## $U(\\theta)$ for several $W_{D3}$"),
        ("code", dedent("""\
            fig, ax = plt.subplots(figsize=(7, 4))
            for W in [0.0, 0.1, 0.5, 1.0]:
                U = full.steering_drive(theta, p.replace(W_D3=W))
                ax.plot(theta, U, label=fr'$W_{{D3}}={W}$')
            ax.axhline(0, color='0.6', lw=0.5)
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$U$')
            ax.legend(); plt.show()
        """)),
        ("md", "## PFL2 population activity peaks at the anti-goal"),
        ("code", dedent("""\
            cont = ContinuousPFLModel(n_quad=256)
            G_P2 = cont.pfl2_population_sum(theta, p)
            fig, ax = plt.subplots(figsize=(7, 3.5))
            ax.plot(theta, G_P2)
            ax.axvline(p.delta, color='C2', ls='--', label=r'$+\\delta$')
            ax.axvline(-p.delta, color='C3', ls='--', label=r'$-\\delta$')
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$\\Sigma_{P2}(\\theta)$')
            ax.legend(); plt.show()
        """)),
    ]


def nb_10():
    return [
        ("md", "# 10 — Extensions\n\n"
               "Asymmetric goal amplitudes, three goals, and parametric "
               "$\\delta$ forcing. (Noise / Kramers escape and slow $\\delta$ "
               "dynamics are sketched but not fully developed.)"),
        ("code", PREAMBLE),
        ("md", "## Asymmetric amplitudes\n\n"
               "Build a thin override of the continuous model with "
               "$A_1 \\neq A_2$ on the two goal bumps. The pitchfork "
               "unfolds into an imperfect bifurcation."),
        ("code", dedent("""\
            from steering.nonlinearities import get_nonlinearity
            from scipy.integrate import fixed_quad

            class AsymmetricContinuous(ContinuousPFLModel):
                def __init__(self, A1=1.0, A2=1.0, n_quad=256):
                    super().__init__(n_quad=n_quad)
                    self.A1 = A1; self.A2 = A2

                def _population_sum(self, theta, params, shift):
                    f = get_nonlinearity(params.nonlinearity, params.nonlinearity_params)
                    theta_arr = np.atleast_1d(np.asarray(theta, dtype=float))
                    phi = self._phi_nodes[None, :]; theta_b = theta_arr[:, None]
                    h = (params.S * np.exp(params.kappa_h * np.cos(theta_b - phi + shift))
                         + self.A1 * np.exp(params.kappa_g * np.cos(params.delta - phi))
                         + self.A2 * np.exp(params.kappa_g * np.cos(-params.delta - phi)))
                    fh = f(h)
                    return np.einsum('ij,j->i', fh, self._phi_weights)

            theta = np.linspace(-np.pi, np.pi, 401)
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            fig, ax = plt.subplots(figsize=(7, 4))
            for r in [1.0, 0.9, 0.7, 0.5]:
                m = AsymmetricContinuous(A1=1.0, A2=r)
                U = m.steering_drive(theta, p)
                ax.plot(theta, U, label=fr'$A_2/A_1={r}$')
            ax.axhline(0, color='0.6', lw=0.5)
            ax.set_xlabel(r'$\\theta$'); ax.set_ylabel(r'$U(\\theta)$')
            ax.legend(); plt.show()
        """)),
        ("md", "## Parametric $\\delta(t)$ forcing"),
        ("code", dedent("""\
            bessel = BesselSteeringModel()
            p = ModelParams(kappa_h=2.0, kappa_g=2.0, delta=1.4)
            dyn = AccelerationDynamics(model=bessel, gamma=0.1)
            forc = ForcingParams(F=0.3, omega=2.0, forcing_type='parametric_delta')
            sim = Simulation(dyn, p, forc, rtol=1e-9, atol=1e-11)
            res = sim.run(np.array([0.5, 0.0]), (0.0, 80.0))
            fig, ax = plt.subplots(figsize=(8, 3.5))
            ax.plot(res.t, res.states[:, 0])
            ax.set_xlabel('t'); ax.set_ylabel(r'$\\theta$')
            ax.set_title(r'Parametric $\\delta(t)$ forcing')
            plt.show()
        """)),
    ]


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------


NOTEBOOKS = {
    "00_model_walkthrough.ipynb": nb_00,
    "01_steering_profiles.ipynb": nb_01,
    "02_pitchfork_bifurcation.ipynb": nb_02,
    "03_phase_portraits_hamiltonian.ipynb": nb_03,
    "04_homoclinic_orbits.ipynb": nb_04,
    "05_melnikov_analysis.ipynb": nb_05,
    "06_poincare_chaos.ipynb": nb_06,
    "07_topology_comparison.ipynb": nb_07,
    "08_discrete_neurons.ipynb": nb_08,
    "09_pfl2_indirect_pathway.ipynb": nb_09,
    "10_extensions.ipynb": nb_10,
}


def build(out_dir: Path | None = None) -> list[Path]:
    out_dir = Path(out_dir) if out_dir is not None else Path(__file__).resolve().parent
    out_dir.mkdir(exist_ok=True, parents=True)
    written: list[Path] = []
    for name, factory in NOTEBOOKS.items():
        nb = nbf.v4.new_notebook()
        cells = []
        for kind, src in factory():
            if kind == "md":
                cells.append(nbf.v4.new_markdown_cell(src))
            elif kind == "code":
                cells.append(nbf.v4.new_code_cell(src))
            else:
                raise ValueError(f"unknown cell kind: {kind}")
        nb["cells"] = cells
        nb["metadata"] = {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python"},
        }
        path = out_dir / name
        with path.open("w") as fh:
            nbf.write(nb, fh)
        written.append(path)
    return written


if __name__ == "__main__":
    paths = build()
    for p in paths:
        print(p)
