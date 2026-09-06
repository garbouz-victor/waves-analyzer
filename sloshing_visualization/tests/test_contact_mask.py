import matplotlib.pyplot as plt
from sloshing.visualization.panels import contact_mask
from sloshing.visualization.scope import load_scope,CONTACT_WARNING


def test_contact_hatching_and_warning_always_added():
    fig,ax=plt.subplots()
    artists=contact_mask(ax,load_scope())
    hatched=[a for a in artists if a.get_hatch()]
    assert len(hatched)==2
    assert any(CONTACT_WARNING==t.get_text() for t in ax.texts)
    spans=sorted((a.get_x(),a.get_x()+a.get_width()) for a in hatched)
    assert spans==[(-1.,-.98),(.98,1.)]
    plt.close(fig)
