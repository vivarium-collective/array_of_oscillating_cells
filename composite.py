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
    return {
        'processes': processes,
        'topology': topology
    }


def run_composite(
        gridr=49,
        gridc=49,
        boundary_molecules=None,
        total_time=10,
):
    # make the composite
    tissue_composite = make_composite(gridr=gridr,
                                      gridc=gridc,
                                      boundary_molecules=boundary_molecules)

    # make the simulation
    tissue = Engine(**tissue_composite)

    # run the simulation
    tissue.update(total_time)

    data = tissue.emitter.get_data()
    print(data)

    # plot
    time_slices = [total_time]
    mol_ids = [
        'X',
        'Y',
        'Yex',
        'Xex'
    ]
    plot_heatmaps(data, time_slices, mol_ids,
                  # output_dir='output', filename='composite'
                  )



if __name__ == '__main__':
    run_composite(total_time=10,
                  gridr=9,
                  gridc=9
                  )
