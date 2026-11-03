#!/usr/bin/env python3
"""Compile the actual AUSM pressure-correction expression with scalar stand-ins.

No OpenFOAM installation is needed. This tests the branch itself, not a second
implementation of the intended formula. Full solver validation is separate.
"""
from pathlib import Path
import subprocess
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[2] / 'src/finiteVolume/fluxSchemes/ausmPlusUp/ausmPlusUpFluxScheme.C'


class PressureCorrectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = SOURCE.read_text()
        start = source.index('\n    Mach_1_2.ref() +=', source.index('// Low Mach number diffusive term'))
        stop = source.index('    MDiff.clear();', start)
        expression = source[start:stop]
        cls.directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.directory.cleanup)
        root = Path(cls.directory.name)
        cpp = root/'probe.cpp'
        cpp.write_text('''#include <iostream>
#include <iomanip>
using scalar = double;
struct Field { double value; double& ref() { return value; } double operator()() const { return value; } };
template<class F> double calcForEachFace(F f, double a, double b) { return f(a,b); }
double corrected(double mach, double pressure) {
    Field Mach_1_2{mach}, MDiff{pressure};
''' + expression + '''
    return Mach_1_2();
}
int main() { double m,d; std::cout << std::setprecision(17);
while(std::cin >> m >> d) std::cout << corrected(m,d) << "\\n"; }
''')
        cls.executable = root/'probe'
        subprocess.run(['c++','-std=c++17',str(cpp),'-o',str(cls.executable)],check=True)

    def evaluate(self, pairs):
        result = subprocess.run([str(self.executable)],input=''.join(f'{m:.17g} {d:.17g}\n' for m,d in pairs),text=True,capture_output=True,check=True)
        return [float(s) for s in result.stdout.split()]

    def test_continuous_when_pressure_cancels_advection(self):
        eps=1e-9
        for sign in (-1,1):
            a,b = self.evaluate([(sign*.1-eps,-sign*.1),(sign*.1+eps,-sign*.1)])
            self.assertLess(abs(b-a),3*eps, f'finite jump at mass-flux reversal: {a}, {b}')

    def test_continuous_when_uncorrected_mach_crosses_zero(self):
        eps=1e-9
        a,b = self.evaluate([(-eps,.1),(eps,.1)])
        self.assertLess(abs(b-a),3*eps)

    def test_full_pressure_term_in_both_flow_directions(self):
        pairs=[(.1,-.1),(-.1,.1),(.1,-.2),(-.1,.2),(.2,.03),(-.2,-.03),(0.,.1)]
        for (m,d),value in zip(pairs,self.evaluate(pairs)):
            self.assertAlmostEqual(value,m+d,places=14)

    def test_zero_pressure_jump_leaves_mach_unchanged(self):
        pairs=[(m,0.) for m in (-2.,-.1,0.,.1,2.)]
        self.assertEqual(self.evaluate(pairs),[m for m,_ in pairs])


if __name__ == '__main__':
    unittest.main()
