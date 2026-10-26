"""Exact analytical examples, not experimental observations or a fitted model.

Run with Python 3.12+. Standard-library numerical analysis; reportlab is used
only for the vector figure. The data contain no measurements from patients.
"""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def endpoint(initial, live, low_per_cell, high_per_cell):
    if initial <= 0 or live < 0:
        raise ValueError('Initial dose must be positive and live count nonnegative')
    low = live * low_per_cell
    high = live * high_per_cell
    return {'initial_beta_cells': initial, 'live_beta_cells': live,
            'basal_output': low, 'stimulated_output': high,
            'incremental_output': high-low, 'retained_abundance_fraction': live/initial,
            'output_per_initial_cell': (high-low)/initial,
            'increment_per_survivor': high_per_cell-low_per_cell if live else None,
            'stimulation_index': high/low if low else None}

def calculate():
    examples = {
        'challenged_control': endpoint(100,40,.5,1.5),
        'persistence_with_lower_output': endpoint(100,80,.2,.6),
        'persistence_with_preserved_per_cell_output': endpoint(100,80,.5,1.5),
        'complete_loss': endpoint(100,0,.5,1.5),
        'secretory_injury_control': endpoint(100,40,.25,.75),
        'secretory_preservation_same_abundance': endpoint(100,40,.5,1.5),
    }
    control = examples['challenged_control']
    for key,e in examples.items():
        reference = 'secretory_injury_control' if key.startswith('secretory_') else 'challenged_control'
        e['reference_scenario'] = reference
        e['output_ratio_to_challenged_control'] = e['output_per_initial_cell']/examples[reference]['output_per_initial_cell']
    checks = {}
    checks['counterexample_same_index'] = math.isclose(examples['persistence_with_lower_output']['stimulation_index'], control['stimulation_index'])
    checks['double_retained_abundance'] = math.isclose(examples['persistence_with_lower_output']['retained_abundance_fraction']/control['retained_abundance_fraction'], 2)
    checks['twenty_percent_less_output'] = math.isclose(examples['persistence_with_lower_output']['output_ratio_to_challenged_control'], .8)
    checks['preserved_per_cell_double_output'] = math.isclose(examples['persistence_with_preserved_per_cell_output']['output_ratio_to_challenged_control'], 2)
    checks['complete_loss_not_dropped'] = examples['complete_loss']['output_per_initial_cell']==0 and examples['complete_loss']['increment_per_survivor'] is None
    checks['factorization_and_frontier'] = all(math.isclose(s*f, (s*control['retained_abundance_fraction']*f*control['increment_per_survivor'])/control['output_per_initial_cell']) and math.isclose(s*(1/s),1) for s in [.25,.5,1,2,4] for f in [.1,.4,1,2,4])
    sc=examples['secretory_injury_control'];se=examples['secretory_preservation_same_abundance']
    checks['function_gain_without_extra_retention'] = sc['live_beta_cells']==se['live_beta_cells']==40 and math.isclose(sc['increment_per_survivor'],.5) and math.isclose(se['increment_per_survivor'],1) and math.isclose(se['output_ratio_to_challenged_control'],2)
    unchallenged_control=endpoint(100,100,.5,1.5)
    unchallenged_edit=endpoint(100,100,.5,1.5)
    interaction=(se['output_per_initial_cell']-sc['output_per_initial_cell'])-(unchallenged_edit['output_per_initial_cell']-unchallenged_control['output_per_initial_cell'])
    checks['matched_unchallenged_function_contrast'] = math.isclose(interaction,.2) and math.isclose(unchallenged_edit['output_per_initial_cell']-unchallenged_control['output_per_initial_cell'],0)
    assert all(checks.values()), checks
    result = {'data_type':'Exact hypothetical examples in arbitrary secretion units; no biological measurements',
              'identity':'F = (QH - QL) / N0 = (Nlive / N0) * ((QH - QL) / Nlive)',
              'domain':'Matched-duration assays, beta-cell-attributable secretion, N0 > 0; survivor factor requires Nlive > 0. Ratios require positive control increments.',
              'examples':examples,
              'secretory_example_unchallenged_controls':{'control':unchallenged_control,'edited':unchallenged_edit,'edit_by_challenge_difference_in_F':interaction},
              'interpretation':'Product output, abundance and survivor-preservation claims are separate. No adequate-output threshold is assigned to these examples.',
              'checks':checks}
    (ROOT/'analysis_results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    with (ROOT/'analytical_examples.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['scenario']+list(control));w.writeheader()
        for key,e in examples.items():w.writerow({'scenario':key,**e})
    return result

def figure(result):
    from reportlab.pdfgen import canvas
    from reportlab.lib.colors import HexColor, Color
    path=ROOT/'figures';path.mkdir(exist_ok=True)
    c=canvas.Canvas(str(path/'survival_function.pdf'),pagesize=(504,320),invariant=1)
    c.setTitle('Retained abundance and per-survivor output ratios')
    c.setAuthor('UCP Technology LLC')
    x0,y0,w,h=62,54,300,230
    def xy(s,f):return x0+(math.log2(s)+2)*w/4,y0+(math.log2(f)+2)*h/4
    c.setFillColor(HexColor('#F1F6F7'));c.rect(x0,y0,w,h,fill=1,stroke=0)
    c.setFillColor(HexColor('#FCF2EE'));p=c.beginPath();p.moveTo(x0,y0);p.lineTo(x0+w,y0);p.lineTo(x0,y0+h);p.close();c.drawPath(p,fill=1,stroke=0)
    c.setFont('Helvetica',9)
    for v in [.25,.5,1,2,4]:
        x,y=xy(v,v);c.setStrokeColor(HexColor('#D6DCE0'));c.setLineWidth(.5)
        c.line(x,y0,x,y0+h);c.line(x0,y,x0+w,y)
        c.setFillColor(HexColor('#333333'));c.drawCentredString(x,y0-14,str(v));c.drawRightString(x0-9,y-3,str(v))
    c.setStrokeColor(HexColor('#283F4B'));c.setLineWidth(1.7);c.line(x0,y0+h,x0+w,y0)
    c.setFillColor(HexColor('#283F4B'));c.setFont('Helvetica-Bold',10)
    c.drawCentredString(252,296,'Output ratio = retained-abundance ratio × per-survivor output ratio')
    c.setFont('Helvetica',10);c.drawCentredString(x0+w/2,17,'Retained viable beta-cell abundance ratio, edited/control')
    c.saveState();c.translate(18,y0+h/2);c.rotate(90);c.drawCentredString(0,0,'Incremental secretion per survivor ratio');c.restoreState()
    c.setFont('Helvetica',9);c.setFillColor(HexColor('#335A60'));c.drawString(236,259,'Higher total output')
    c.setFillColor(HexColor('#8C4B36'));c.drawString(77,72,'Lower total output')
    points=[('A',1,1,'Reference','#333333'),('B',2,.4,'0.8 x output','#B94C35'),('C',2,1,'2 x output','#147A78'),('D',1,2,'2 x output','#72549B')]
    for label,s,f,desc,color in points:
        x,y=xy(s,f);c.setFillColor(HexColor(color));c.circle(x,y,4,fill=1,stroke=0);c.setFont('Helvetica-Bold',10);c.drawString(x+7,y+4,label)
    c.setFont('Helvetica-Bold',9);c.setFillColor(HexColor('#222222'));c.drawString(379,263,'Illustrative points')
    c.setFont('Helvetica',9)
    for i,(label,s,f,desc,color) in enumerate(points):
        yy=242-i*44;c.setFillColor(HexColor(color));c.drawString(379,yy,label+'  '+desc)
        c.setFillColor(HexColor('#333333'));c.drawString(379,yy-13,'Abundance '+str(s)+' x')
        c.drawString(379,yy-26,'Per cell '+str(f)+' x')
    c.setFont('Helvetica',8);c.drawString(379,53,'Line: equal total output')
    c.drawString(379,40,'No observed data plotted')
    c.save()

if __name__=='__main__':
    r=calculate();figure(r)
    print(json.dumps({'analytical_checks_passed':sum(r['checks'].values()),'empirical_observations':0}))
