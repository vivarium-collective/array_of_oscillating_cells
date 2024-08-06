import random
import time

from process import CellProcess, DEFAULT_MODEL_FILE
from plotting import plot_heatmaps
from vivarium.core.engine import Engine
from basico import load_model


def make_composite(
    gridr=49,
    gridc=49,
    boundary_molecules=None,
):
    if boundary_molecules is None:
        boundary_molecules = ['Xex']

    # Load the single cell model into Basico
    copasi_model_object = load_model(DEFAULT_MODEL_FILE)

    # make the composite
    processes = {}
    topology = {}
    initial_state = {}
    for row in range(0, gridr):
        for col in range(0, gridc):
            config = {
                'copasi_object': copasi_model_object,
                'boundary_molecules': boundary_molecules,
                'n_sides': 4,
                'time_step': 1
            }
            cell_id = f'[{row},{col}]'
            cell = CellProcess(config)
            processes[cell_id] = cell
            topology[cell_id] = {
                'boundary': (f'{cell_id}_store', 'boundary'),
                'internal': (f'{cell_id}_store', 'internal')
            }

            initial_state[f'{cell_id}_store'] = {
                'boundary': {
                    mol_id: random.uniform(0, 0.1) for mol_id in boundary_molecules
                },
            }
    return {
        'processes': processes,
        'topology': topology,
        'initial_state': initial_state,
    }


def run_composite(
        gridr=49,
        gridc=49,
        boundary_molecules=None,
        total_time=10,
):
    # make the composite

    # time how long it takes
    start = time.time()
    tissue_composite = make_composite(gridr=gridr,
                                      gridc=gridc,
                                      boundary_molecules=boundary_molecules)
    print(f"Time to initialize composite: {time.time() - start}")

    # make the simulation
    tissue = Engine(**tissue_composite, display_info=True)

    # run the simulation
    start = time.time()
    tissue.update(total_time)
    print(f"Time to run simulation: {time.time() - start}")

    data = tissue.emitter.get_data()
    # print(data)

    # plot
    time_slices = [total_time]
    mol_ids = [
        'X',
        'Y',
        'Yex',
        'Xex'
    ]
    plot_heatmaps(data, time_slices, mol_ids,
                  output_dir='output', filename='composite'
                  )



if __name__ == '__main__':
    run_composite(total_time=10,
                  gridr=49,
                  gridc=49
                  )
