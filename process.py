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

def _set_parameters(changes, dm):
    model = dm.getModel()
    assert (isinstance(model, COPASI.CModel))

    references = COPASI.ObjectStdVector()

    for name, value in changes:
        parameter = model.getModelParameterSets().getParameter(name)
        assert (isinstance(parameter, COPASI.CModelValue))
        if parameter is None:
            print(f"Parameter {name} not found in model")
            continue
        parameter.setValue(value)
        references.append(parameter.getValueReference())

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
        'parameter_noise': {'alpha': 0.1},
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

        # get the parameters for which we need to add noise
        self.noise_parameters = {}
        for param, noise in self.parameters['parameter_noise'].items():
            p = get_parameters(name=param, model=self.copasi_model_object)
            value = p['initial_value'][0]
            # add noise to the parameter from uniform distribution
            sampled = np.random.uniform(-noise, noise)
            self.noise_parameters[param] = value + sampled

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

        # set species concentrations
        changes = []
        for mol_id, value in states['boundary'].items():
            # if mol_id.endswith('_ext'):
            changes.append((mol_id, value))
        for mol_id, value in states['internal'].items():
            changes.append((mol_id, value))
        _set_initial_concentrations(changes, self.copasi_model_object)

        # set parameters
        for param, value in self.noise_parameters.items():
            set_reaction_parameters(name=param, value=value, model=self.copasi_model_object)
            # self.copasi_model_object.set_reaction_parameters(param, value)

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
