"""Heggie-only plotting, metadata, sampling and CLI checks."""
import sys
from pathlib import Path
import tempfile
import subprocess
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"script"))
from heggie_reference import sample
from plot_heggie_hydro import load_run, load_stock_run, segregation_curves, diagnostics
from plot_heggie_calibrated import density_table


class HeggieTest(unittest.TestCase):
    def test_sampling(self):
        data=np.ones((3,13));data[:,1]=[0,1,2];data[:,2]=[1,2,3]
        np.testing.assert_allclose(sample(data,[.5]),[[1.5,1]])
        with self.assertRaisesRegex(ValueError,"no extrapolation"):sample(data,[3])

    def test_table(self):
        text=density_table([[1,3.6,.6]],[[1.35,.4]],[[3.6,.6]])
        self.assertTrue(text.startswith(r"\begin{tabular}"))
        self.assertTrue(text.endswith("\\bottomrule\n\\end{tabular}\n"))

    def test_growth_rate_units(self):
        data=np.ones((101,13));data[:,0]=np.linspace(0,1,101)
        data[:,2]=np.exp(2*data[:,0]);data[:,3]=np.exp(3*data[:,0])
        data[:,9]=4;data[:,10]=5
        _,_,x,y=segregation_curves(data)
        np.testing.assert_allclose(x,8,rtol=1e-12)
        np.testing.assert_allclose(y,15,rtol=1e-12)
        self.assertNotIn("max_step_energy_ledger_error",diagnostics(data,"TEST"))

    def test_history_validation(self):
        (ROOT/"output").mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT/"output") as tmp:
            root=Path(tmp);data=np.ones((3,13));data[:,0]=[0,1,2]
            data.astype("<f8").tofile(root/"history.dat")
            (root/"status.txt").write_text("TIME_LIMIT\n")
            np.testing.assert_array_equal(load_run(root)[0],data)
            (root/"solver.txt").write_text("hydrostatic\n")
            with self.assertRaisesRegex(ValueError,"unchanged-solver"):load_stock_run(root)
            (root/"solver.txt").write_text("hydrostatic-stock\n")
            np.testing.assert_array_equal(load_stock_run(root)[0],data)
            data[1,0]=0;data.astype("<f8").tofile(root/"history.dat")
            with self.assertRaisesRegex(ValueError,"unordered"):load_run(root)

    @unittest.skipUnless((ROOT/"main-strict").exists(),"build main-strict for CLI/readback integration")
    def test_cli_and_serialization(self):
        import os
        with tempfile.TemporaryDirectory(dir=ROOT/"output") as tmp:
            out=Path(tmp)/"run"
            command=[str(ROOT/"main-strict"),"hydro-heggie",str(out),"3",".005","60","1e-8",".001","20",".72","1.5"]
            env={**os.environ,"OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1"}
            subprocess.run(command,check=True,capture_output=True,env=env)
            d,status=load_run(out);self.assertEqual(status,"TIME_LIMIT")
            self.assertEqual((out/"solver.txt").read_text().splitlines()[0],"hydrostatic-stock")
            np.testing.assert_array_equal(np.fromfile(out/"profile/param.dat",dtype="<f8"),[.72,1.5])
            self.assertEqual(np.fromfile(out/"diagnostics.dat",dtype="<f8").size,8*len(d))
            original=(out/"history.dat").read_bytes()
            self.assertNotEqual(subprocess.run(command,capture_output=True,env=env).returncode,0)
            self.assertEqual((out/"history.dat").read_bytes(),original)
