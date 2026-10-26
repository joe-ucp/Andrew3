"""Reproduce the analytical example, figure and Word manuscript.

Requires Python 3.12+, the declared packages, and Poppler pdftoppm on PATH.
PDF manuscript export is a separate documented LibreOffice command.
"""
from pathlib import Path
import json
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
    for script in ['analysis.py','write_content.py']:
        subprocess.run([sys.executable,str(ROOT/script)],check=True,cwd=ROOT)
    poppler=shutil.which('pdftoppm')
    if not poppler:
        raise SystemExit('Poppler pdftoppm is required to reproduce the figure PNG. No PDF export was performed.')
    subprocess.run([poppler,'-singlefile','-r','180','-png',str(ROOT/'figures/survival_function.pdf'),str(ROOT/'figures/survival_function')],check=True)
    subprocess.run([sys.executable,str(ROOT/'build_manuscript.py')],check=True,cwd=ROOT)
    results=json.loads((ROOT/'analysis_results.json').read_text(encoding='utf-8'))
    print(json.dumps({'status':'PASS','analytical_checks':sum(results['checks'].values()),'biological_experiments':0,'PDF_manuscript_compiled':False}))
