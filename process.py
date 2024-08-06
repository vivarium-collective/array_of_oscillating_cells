"""
https://github.com/copasi/model_replicator/tree/main/examples/Array_of_oscillating_cells
"""
from vivarium.core.process import Process
from vivarium.core.engine import Engine, pp
from basico import *

import COPASI

DEFAULT_MODEL_FILE = 'glycolysis-autocatalytic.cps'


def _set_initial_concentrations(changes, dm):
    model = dm.getModel()
    assert (isinstance(model, COPASI.CModel))

    references = COPASI.ObjectStdVector()

    for name, value in changes:
        species = model.getMetabolite(name)
        assert (isinstance(species, COPASI.CMetab))
        if species is None:
            print(f"Species {name} not found in model")
            continue
        species.setInitialConcentration(value)
        references.append(species.getInitialConcentrationReference())

    model.updateInitialValues(references)


def _get_transient_concentration(name, dm):
    model = dm.getModel()
    assert (isinstance(model, COPASI.CModel))
    species = model.getMetabolite(name)
    assert (isinstance(species, COPASI.CMetab))
    if species is None:
        print(f"Species {name} not found in model")
        return None
    return species.getConcentration()


class CellProcess(Process):
    defaults = {
        'model_file': DEFAULT_MODEL_FILE,
        'copasi_object': None,
        'boundary_molecules': ['Xex'],
        'n_sides': 4,
        'time_step': 1.0,
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

        # Load the single cell model into Basico
        if self.parameters['copasi_object']:
            self.copasi_model_object = self.parameters['copasi_object']
        else:
            self.copasi_model_object = load_model(self.parameters['model_file'])
        all_species = get_species(model=self.copasi_model_object).index.tolist()
        self.internal_species = [
            species for species in all_species if species not in self.parameters['boundary_molecules']]
        self.boundary_species = self.parameters['boundary_molecules']

    def ports_schema(self):
        ports = {
            'boundary': {
                mol_id: {
                    '_default': _get_transient_concentration(name=mol_id, dm=self.copasi_model_object),
                    '_emit': True,
                    '_updater': 'set',
                } for mol_id in self.boundary_species
            },
            'internal': {
                mol_id: {
                    '_default': _get_transient_concentration(name=mol_id, dm=self.copasi_model_object),
                    '_emit': True,
                    '_updater': 'set',
                } for mol_id in self.internal_species
            }
        }

        return ports

    def next_update(self, endtime, states):

        # set boundary species concentrations
        boundary_species = states['boundary']
        changes = []
        for mol_id, value in boundary_species.items():
            # if mol_id.endswith('_ext'):
            changes.append((mol_id, value))

        _set_initial_concentrations(changes, self.copasi_model_object)

        # run model for "endtime" length; we only want the state at the end of endtime, if we need more we can set intervals to a larger value
        timecourse = run_time_course(duration=endtime, intervals=1, update_model=True, model=self.copasi_model_object)

        # extract end values of concentrations from the model and set them in results (18 states)
        results = {'boundary': {}, 'internal': {}}
        for mol_id in self.boundary_species:
            results['boundary'][mol_id] = _get_transient_concentration(name=mol_id, dm=self.copasi_model_object)
        for mol_id in self.internal_species:
            results['internal'][mol_id] = _get_transient_concentration(name=mol_id, dm=self.copasi_model_object)

        return results


def test_single_cell():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    model_path = os.path.join(current_dir, f"{DEFAULT_MODEL_FILE}")
    total_time = 1000
    config = {
        'model_file': model_path,
        'boundary_molecules': ['Xex'],
        'n_sides': 4,
        'time_step': 100
    }
    cell = CellProcess(config)

    ports = cell.ports_schema()
    print('PORTS')
    print(ports)

    sim = Engine(
        processes={'cell': cell},
        topology={'cell': {port_id: (port_id,) for port_id in ports.keys()}}
    )

    sim.update(total_time)

    data = sim.emitter.get_timeseries()
    print('RESULTS')
    print(pp(data))


if __name__ == '__main__':
    test_single_cell()
