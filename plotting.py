import matplotlib.pyplot as plt
import numpy as np
import os


def plot_heatmaps(data, time_slices, mol_ids, output_dir=None, filename=None):
    for time_index in time_slices:
        if time_index not in data:
            continue

        time_data = data[time_index]

        # Extract the row and column size from the cell names
        cells = list(time_data.keys())
        max_row = max(int(cell.split(',')[0][1:]) for cell in cells) + 1
        max_col = max(int(cell.split(',')[1].split(']')[0]) for cell in cells) + 1

        for mol_id in mol_ids:
            heatmap_internal = np.zeros((max_row, max_col))
            heatmap_boundary = np.zeros((max_row, max_col))
            found_internal = False
            found_boundary = False

            for cell, cell_data in time_data.items():
                row = int(cell.split(',')[0][1:])
                col = int(cell.split(',')[1].split(']')[0])

                if mol_id in cell_data['internal']:
                    heatmap_internal[row, col] = cell_data['internal'][mol_id]
                    found_internal = True

                if mol_id in cell_data['boundary']:
                    heatmap_boundary[row, col] = cell_data['boundary'][mol_id]
                    found_boundary = True

            if found_internal:
                plt.figure(figsize=(8, 6))
                plt.imshow(heatmap_internal, cmap='viridis', interpolation='none')
                plt.colorbar()
                plt.title(f'Internal {mol_id} at Time Index {time_index}')
                plt.xlabel('Column')
                plt.ylabel('Row')

                if output_dir and filename:
                    file_path = os.path.join(output_dir, f"{filename}_internal_{mol_id}_time_{time_index}.png")
                    plt.savefig(file_path)
                else:
                    plt.show()

            if found_boundary:
                plt.figure(figsize=(8, 6))
                plt.imshow(heatmap_boundary, cmap='viridis', interpolation='none')
                plt.colorbar()
                plt.title(f'Boundary {mol_id} at Time Index {time_index}')
                plt.xlabel('Column')
                plt.ylabel('Row')

                if output_dir and filename:
                    file_path = os.path.join(output_dir, f"{filename}_boundary_{mol_id}_time_{time_index}.png")
                    plt.savefig(file_path)
                else:
                    plt.show()

