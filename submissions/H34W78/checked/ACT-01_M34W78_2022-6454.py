import math as mt
import scipy.optimize as spo

#Load 1: 400 kW at 0.85 pf lagging
P1 = 400
pf1 = 0.85
Q1 = P1 * mt.tan(mt.acos(pf1))

#Load 2: 300 kVAR at 0.75 pf lagging
Q2 = 300
pf2 = 0.75
P2 = Q2 / mt.tan(mt.acos(pf2))

S_total = complex(P1 + P2, Q1 + Q2) * 1000  #[VA]

Z = complex(1.0, 2.5)

def resid_from_x(V):
    V_rms = V * 1000  #convert kV -> V
    I = S_total.conjugate() / V_rms
    V_send = V_rms + I * Z
    return abs(V_send) / 1000 - 13.8  #result in kV

def find_valid_bounds(f, start=1, stop=50, step=0.5):
    a = start
    while a < stop:
        b = a + step
        if f(a) * f(b) < 0:
            return a, b
        a += step
    raise ValueError("No valid bounds found where function changes sign.")

XL, XU = find_valid_bounds(resid_from_x)
print(f"✅ Valid interval found: [{XL:.2f}, {XU:.2f}] kV")

X_TOL = 1e-7
MAX_ITERS = 750

p_bs, info = spo.bisect(
    resid_from_x,
    XL,
    XU,
    xtol=X_TOL,
    maxiter=MAX_ITERS,
    full_output=True
)

print("\n--- BISECTION METHOD RESULTS ---")
print(f"Receiving-end Voltage (Root) : {p_bs:.6f} kV")
print(f"Iterations                   : {info.iterations}")
print(f"Converged?                   : {info.converged}")
print(f"Final Residual (kV)          : {resid_from_x(p_bs):.6e}")
