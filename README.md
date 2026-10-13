# Chimera States and Plasticity

Simulation code accompanying the peer-reviewed publication:

> Yang, C.; Santos, M.S.; **Protachevicz, P.R.**; Reis, P.D.C.; et al. **"Chimera states induced by spike timing-dependent plasticity in a regular neuronal network."** *AIP Advances*, 12, 105119, 2022. DOI: [10.1063/5.0101055](https://doi.org/10.1063/5.0101055)

Originally published as part of the author's research data collection at [henon.if.usp.br/OscilControlData](http://henon.if.usp.br/OscilControlData/).

## Research summary

Chimera states are spatiotemporal patterns in which distinct dynamics coexist, such as synchronous and asynchronous firing in the same network. This work studies how spike timing-dependent plasticity (STDP) affects the emergence and lifetime of chimera states in a **regular network of coupled adaptive exponential integrate-and-fire neurons**, where all connections initially carry the same synaptic strength, which STDP then reshapes based on the timing between pre- and postsynaptic spikes.

Key findings:
- The parameter range producing chimera states is **larger** in the plastic network than in a network without plasticity
- Chimera **lifetime increases** when STDP is active
- Neuronal spike frequency increases under constant positive current
- Changes in synaptic weights induced by STDP increase the overall appearance of chimera states in parameter space

## What the code does

Two complementary simulation programs are provided, both built on the same adaptive exponential integrate-and-fire (AdEx) neuron model with STDP-governed synaptic weights:

- **`plastic110.c`**: sweeps the external current `I0` and the coupling range `r` to map out where chimera, synchronized, and desynchronized regimes occur, tracking the mean synaptic weight (`gexM`), coefficient of variation, and firing frequency across the parameter space. Reads an index from stdin to parallelize runs across `I0` ranges (`I0_r_<variable>.dat` output).
- **`plastic_potentiation_matrix.c`**: focuses on the evolution of the full synaptic weight matrix for specific initial conditions (fake-chimera, tilted, and synchronized patterns), saving the raster of spikes, the final weight matrix (`pesos_*.dat`), and summary potentiation statistics (`diags_*.dat`) per run. `compilar_rodar.sh` compiles and launches 4 parallel runs (one per initial-condition index).

Both share the same core structure: a ring/band-coupled network of `N≈1000` neurons, STDP-governed synaptic updates (potentiation/depression based on spike-timing deltas), and 4th-order Runge-Kutta integration of the AdEx membrane dynamics.

## Requirements

A C compiler (GCC, ICC, or Intel's `icx`) and `libm`. No external dependencies.

## Build & run

```bash
# Parameter-space sweep
gcc plastic110.c -lm -O2 -o chimera_sweep
echo 1 | ./chimera_sweep

# Weight-matrix / initial-condition study (parallel runs)
bash compilar_rodar.sh
```

Output `.dat` files can be plotted with Gnuplot, Python/Matplotlib, or any standard plotting tool.
