---
name: debruijn-newman-chladni-evolution
description: Scientific knowledge skill for Paper XI on Riemann's Xi-function under De Bruijn-Newman evolution (tau = t_DBN). Contains the exact Heat-Flow PDE Hierarchy on Jensen-Laguerre invariants (d_tau + 1/8 d_t^2) L_m = ((m+1)(2m+1)/4) L_{m+1}, the corrected (m-1)^2/2 multiple-zero coalescence lift, the nodal-curvature interpretation of the Csordas-Smith-Varga / Rodgers-Tao zero velocity, the boundary kernel zero 2 eta_0 < pi/4, and the conditional pitchfork saddle perturbation.
---

# Heat-Flow PDE Hierarchy on Jensen Invariants, Nodal-Curvature Zero Velocity, and Conditional Saddle Perturbation under De Bruijn–Newman Evolution (Paper XI)

## 1. Scope, Canonical Normalization, and Red-Team Guardrails

This skill encapsulates the verified differential identities, symbolic `SymPy` proofs, and $50$-decimal `mpmath` checks of **Paper XI** (author: Thijs Missiaen, ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)). Always enforce the following canonical normalizations and audit guardrails:

1. **Exact De Bruijn–Newman Time Identification ($\tau = t _ {\mathrm{DBN}}$)**: Comparing $\Xi _ \tau(z) := 8\int _ 0^\infty \Phi(u)\mathrm{e}^{\tau u^2}\cos(2zu) \mathrm{d}u$ with $H _ \tau(x) := \int _ 0^\infty \mathrm{e}^{\tau u^2}\Phi(u)\cos(xu) \mathrm{d}u$ (Rodgers & Tao 2020) shows $\Xi _ \tau(z) = 8H _ \tau(2z)$ with the **identical** integration variable $u$. Hence $\tau = t _ {\mathrm{DBN}}$ (never $4t _ {\mathrm{DBN}}$), $\Lambda _ \tau = \Lambda$, and the prefactor in $\partial _ \tau\Xi _ \tau = -\frac{1}{4}\partial _ z^2\Xi _ \tau$ comes solely from $x = 2z$.
2. **Canonical Wronskian and Non-Zero Jensen Index Limit**: Across all papers, $W _ \tau(d,t) := \frac{1}{2}\partial _ d\lvert\Xi _ \tau(t - \mathrm{i}d)\rvert^2 = 16\int _ 0^\infty y\sinh(2dy) I _ \tau(y,t) \mathrm{d}y$ and $c _ {2m-1}(t;\tau) = m L _ m(\Xi _ \tau; t) = \frac{8}{(2m-1)!}\int _ 0^\infty (2y)^{2m}I _ \tau(y,t) \mathrm{d}y$, where $L _ m(f; t)$ sums over $0 \le j \le 2m$ (an odd upper limit $2m - 1$ vanishes identically by antisymmetry).
3. **Domain & Sign Hypotheses**: At $y = 0$, $G _ {0,\tau}(\mathrm{i}\eta) = \mathrm{e}^{-\frac{\tau}{2}\eta^2}\lvert\Phi(\mathrm{i}\eta/2)\rvert^2$ vanishes at $2\eta _ 0 = 0.6388300536\dots < \frac{\pi}{4}$, so $\ln G _ {y,\tau}(\mathrm{i}\eta)$ requires restricting to zero-free subdomains $\Phi(\frac{y+\mathrm{i}\eta}{2}) \neq 0$, and $\partial _ \tau\eta _ c(\omega; 0) < 0$ is conditional on $1 - (\partial _ y L _ {y _ c}''(\eta _ c))\partial _ \tau y _ c(\omega; 0) > 0$. Backward Lehmer-pair coalescence ($\tau < 0$) belongs to the Csordas–Norfolk–Varga (1988) *lower-bound* regime $\Lambda \ge t _ {\mathrm{lb}}$ and never implies $\Lambda \le 0$.

## 2. Verified Mathematical Identities

### 2.1 Nodal-Curvature Interpretation of the CSV / Rodgers–Tao Velocity (Proposition 2.1)

For simple zeros $\gamma _ k$ ($\Xi(\gamma _ k) = 0$, $\Xi'(\gamma _ k) \neq 0$) and critical-line extrema $\tau _ k$ ($\Xi'(\tau _ k) = 0$, $\Xi''(\tau _ k) \neq 0$), the Csordas–Smith–Varga (1994) and Rodgers–Tao (2020) velocity formulas satisfy

$$
\gamma _ k'(0) = \frac{\Xi''(\gamma _ k)}{4 \Xi'(\gamma _ k)} = \frac{1}{4} \kappa(\gamma _ k) = \frac{1}{2} \mathrm{P.V.} \sum _ {\gamma _ j \neq \gamma _ k} \frac{1}{\gamma _ k - \gamma _ j}, \qquad \tau _ k'(0) = \frac{\Xi'''(\tau _ k)}{4 \Xi''(\tau _ k)} = \frac{3}{4} \kappa _ \Gamma(\tau _ k),
$$

where $\kappa(\gamma _ k) := t _ {\mathcal{R} _ k}''(0)$ and $\kappa _ \Gamma(\tau _ k) := \tau _ {\Gamma _ k}''(0)$ are the parabolic curvatures of $R(d,t) = 0$ and $I(d,t)/d = 0$ at $d = 0$.

### 2.2 Mirror Gaussian Factorization and Conditional Saddle Perturbation (Propositions 3.1 & 3.2)

In mirror coordinates $y = u + v$, $x = u - v$, $(\frac{y+x}{2})^2 + (\frac{y-x}{2})^2 = \frac{1}{2}(y^2 + x^2)$ yields

$$
G _ {y,\tau}(x) = \mathrm{e}^{\frac{\tau}{2}y^2} \mathrm{e}^{\frac{\tau}{2}x^2} G _ {y,0}(x), \qquad \frac{\partial I _ \tau}{\partial\tau}(y,t) = \frac{1}{8}\left(4y^2 - \frac{\partial^2}{\partial t^2}\right)I _ \tau(y,t).
$$

On a zero-free neighborhood $\Phi(\frac{y+\mathrm{i}\eta}{2}) \neq 0$ of a pitchfork saddle $(y _ c(\omega; 0), \eta _ c(\omega; 0))$ ($L _ {y _ c}'(\eta _ c) = \omega$, $L _ {y _ c}''(\eta _ c) = 0$ with $L _ y(\eta) := 2\ln\lvert\Phi(\frac{y+\mathrm{i}\eta}{2})\rvert$):

$$
\frac{\partial y _ c}{\partial\tau}(\omega; 0) = \frac{\eta _ c(\omega; 0)}{\partial _ y L _ {y _ c}'(\eta _ c(\omega; 0))}, \qquad \frac{\partial\eta _ c}{\partial\tau}(\omega; 0) = \frac{1 - (\partial _ y L _ {y _ c}''(\eta _ c)) \partial _ \tau y _ c(\omega; 0)}{L _ {y _ c}'''(\eta _ c(\omega; 0))}.
$$

### 2.3 Second-Order Radial Wave-Heat PDE and Jensen PDE Hierarchy (Theorem 4.1 & Corollary 4.3)

For any real entire family $f(z,\tau)$ solving $\partial _ \tau f = -\frac{1}{4}\partial _ z^2 f$, $M(d,t;\tau) := \lvert f(t - \mathrm{i}d,\tau)\rvert^2 = \sum _ {m=0}^\infty L _ m(f; t) d^{2m}$ and $\widetilde{W} _ f := \frac{1}{2d}\partial _ d M = \sum _ {m=1}^\infty c _ {2m-1}(t;\tau) d^{2m-2}$ ($c _ {2m-1} = m L _ m$) satisfy the second-order PDE $\partial _ \tau\widetilde{W} _ f + \frac{1}{8}(\partial _ t^2 - \partial _ d^2 - \frac{2}{d}\partial _ d)\widetilde{W} _ f = 0$ and the exact hierarchy

$$
\left(\frac{\partial}{\partial\tau} + \frac{1}{8}\frac{\partial^2}{\partial t^2}\right) L _ m(f; t) = \frac{(m + 1)(2m + 1)}{4} L _ {m+1}(f; t) \quad (m \ge 0), \qquad \left(\frac{\partial}{\partial\tau} + \frac{1}{8}\frac{\partial^2}{\partial t^2}\right) c _ {2m-1} = \frac{m(2m + 1)}{4} c _ {2m+1} \quad (m \ge 1).
$$

At an $m$-fold real zero ($m \ge 2$), $\partial _ \tau L _ {m-1}(t _ 0;\tau _ {\mathrm{coal}}) = \frac{m(m-1)}{2}L _ m(t _ 0;\tau _ {\mathrm{coal}}) > 0$ and $\partial _ \tau c _ {2m-3}(t _ 0;\tau _ {\mathrm{coal}}) = \frac{(m-1)^2}{2}c _ {2m-1}(t _ 0;\tau _ {\mathrm{coal}}) > 0$. Across close pairs with half-gap $\delta = \frac{1}{2}(\gamma _ {k+1} - \gamma _ k)$, $\partial _ \tau c _ 1(t _ {\mathrm{mid}}; 0) = \frac{1}{4}(\Xi''(t _ {\mathrm{mid}}))^2 + \frac{1}{4}\Xi(t _ {\mathrm{mid}})\Xi^{(4)}(t _ {\mathrm{mid}})$ and $t _ {\mathrm{DBN,dip}}(t _ {\mathrm{mid}}) = -2c _ 1(t _ {\mathrm{mid}};0)/c _ 3(t _ {\mathrm{mid}};0) = -2\delta^2(1 + \mathcal{O}(\delta^2))$.

## 3. Reproducible Verification Workflow

Run the standalone symbolic and numerical verification suite from the repository root:

```bash
python scripts/verify_paper11_debruijn_newman_chladni.py
```
