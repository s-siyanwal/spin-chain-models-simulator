from spinchain.paper import make,audit

def test_fig2_parameters():
    assert make('polymer').pe==200 and make('weak_spiral').pe==1000 and make('strong_spiral').pe==5000
    assert make('strong_spiral').flexure==25000
def test_tau_and_steps():
    p=make('polymer');assert p.tau==252500;assert p.steps_for_tau(5)==25250000000
def test_audit_honesty():
    a=audit(make('polymer'),5,3,100000);assert a['total_steps']==75750000000 and a['estimated_hours']>100
