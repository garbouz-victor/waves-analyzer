from pathlib import Path
import numpy as np

from sloshing.visualization.explorer import display_numbers
from sloshing.visualization.panels import display_cell_edges


def test_html_quantization_is_pure_and_retains_extrema_not_clipping():
    values=np.array([[-1.23456789e-8,0.,100.1234567]])
    original=values.copy()
    rounded=np.array(display_numbers(values))
    np.testing.assert_array_equal(values,original)
    np.testing.assert_allclose(rounded,values,rtol=5e-6,atol=0)
    assert rounded.min()<0 and rounded.max()>100


def test_html_heatmap_edges_stay_in_reference_domain():
    nodes=np.linspace(-1.5,0,61)
    edges=display_cell_edges(nodes)
    assert len(edges)==62 and edges[-1]==0 and edges[0]==-1.5
    assert np.all(np.diff(edges)>0)


def test_explorer_contract_has_fixed_limits_and_permanent_warning():
    import sloshing.visualization.explorer as module
    template=Path(module.__file__).with_name('explorer_template.html').read_text()
    assert 'zmin:lim[0],zmax:lim[1],zauto:false' in template
    assert 'id="warning"' in template and 'even when shading is switched off' in template
    assert 'M.particle_displacement_magnification*(q-D.seeds[k][j])' in template
    assert 'no interpolated PDE frames' in template
