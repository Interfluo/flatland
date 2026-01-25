import numpy as np
import pyvista as pv


def save_field(filename, values):
    lines = [f"{i}\n" for i in values]
    
    with open(filename, "w") as f:
        f.writelines(lines)
    return None


if __name__ == "__main__":
    mesh = pv.read("bunny.obj")
    
    for i, v in enumerate(np.linspace(4, 20, 21)):
        node_field = np.sin(np.pi*v*mesh.points[:, 0])
        face_field = np.sin(np.pi*v*mesh.cell_centers().points[:, 0])
        
        save_field(f"data/node_field_{i}.txt", node_field)
        save_field(f"data/face_field_{i}.txt", face_field)
