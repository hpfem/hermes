from hermes2d.modules.schroedinger import (ModuleSchroedinger,
        PotentialHarmonicOscillator)

def test_oscillator():
    m = ModuleSchroedinger()
    p = PotentialHarmonicOscillator()
    p.set_omega(5)
    m.set_potential(p)
