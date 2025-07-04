#!/usr/bin/env -S uv run --script --verbose
# /// script
# requires-python = ">=3.12"
# dependencies = [
#    "rasterio"
# ]
# ///

# https://forum.bambulab.com/t/how-do-you-do-the-new-obj-with-mtl-feature/75622/8

from dataclasses import dataclass
import math
from typing import Optional
import random
import rasterio
import cProfile


@dataclass
class Point3:
    x: float
    y: float
    z: float


@dataclass
class Point6:
    x: float
    y: float
    z: float
    r: float
    g: float
    b: float


@dataclass
class Point4:
    x: float
    y: float
    z: float
    w: float


@dataclass
class Face:
    vertex_indices: list[int]
    vertex_normal_indices: list[int]
    texture_vertex_indices: list[int]


@dataclass
class TexCoord:
    u: float
    v: Optional[float] = None
    w: Optional[float] = None


@dataclass
class RGB:
    r: float
    g: float
    b: float


@dataclass
class Material:
    colour: RGB


@dataclass
class Object:
    verticies: list[Point3]
    vertex_normals: list[Point3]
    faces: list[Face]
    texture_verticies: list[TexCoord]
    name: str
    material: Material

# Define function to load geotiff
def read_geotiff(path):
    tiff:rasterio.DatasetReader = rasterio.open(path)
    return tiff

def get_value_for_point(tiff, x, y, band):
    """
    Get the value for a point in the geotiff.
    """
    # Read the pixel value at (x, y)
    return tiff.read(band)[y, x]

def get_colour_for_point(tiff, x, y):
    """
    Get the colour for a point in the geotiff.
    """
    # Read the pixel value at (x, y)
    values  = tiff.sample([(x, y)],[1,2,3])
    return (values[0], values[1], values[2])

def wiggle(i, j):
    return 1.0 * math.sin(i / 20.0) + 1.0 * math.sin(j / 20.0)


def create_surface():
    """
    This function defines a 2d flat surface of 100x100 points, ((2*99*99) faces)
    """
    verticies: list[Point3] = []
    faces: list[Face] = []

    xlen = 1000
    ylen = 1000
    texture = []
    texture.append(TexCoord(0.0, 0.0))
    texture.append(TexCoord(1.0, 0.0))
    texture.append(TexCoord(0.0, 1.0))
    texture.append(TexCoord(1.0, 1.0))

    tiff = read_geotiff("test.tif")
    ele = read_geotiff("test3.tif")

    # Create top surface for the obj
    def create_surface(top: bool, height_offset: float, create_faces=True):
        vertex_offset = len(verticies)
        for i in range(0, xlen):
            for j in range(0, ylen):
                verticies.append(
                    Point6(
                        i,
                        j,
                        height_offset + (get_value_for_point(ele, i, j, 1)/2.0 if top else 0.0),
                        get_value_for_point(tiff,i,j,1)/255.0 if top else 255.0,
                        get_value_for_point(tiff,i,j,2)/255.0 if top else 255.0,
                        get_value_for_point(tiff,i,j,3)/255.0 if top else 255.0,
                    )
                )
                # Join bl, tr, br
                if i < xlen - 1 and j < ylen - 1 and create_faces:
                    if not top:
                        faces.append(
                            Face(
                                vertex_indices=[
                                    vertex_offset + i * ylen + j,
                                    vertex_offset + i * ylen + j + 1,
                                    vertex_offset + (i + 1) * ylen + j,
                                ],
                                vertex_normal_indices=[],
                                texture_vertex_indices=[0, 1, 2],
                            )
                        )
                        faces.append(
                            Face(
                                vertex_indices=[
                                    vertex_offset + (i + 1) * ylen + j,
                                    vertex_offset + i * ylen + j + 1,
                                    vertex_offset + (i + 1) * ylen + j + 1,
                                ],
                                vertex_normal_indices=[],
                                texture_vertex_indices=[2, 1, 3],
                            )
                        )
                    else:
                        faces.append(
                            Face(
                                vertex_indices=[
                                    vertex_offset + i * ylen + j,
                                    vertex_offset + (i + 1) * ylen + j,
                                    vertex_offset + i * ylen + j + 1,
                                ],
                                vertex_normal_indices=[],
                                texture_vertex_indices=[0, 1, 2],
                            )
                        )
                        faces.append(
                            Face(
                                vertex_indices=[
                                    vertex_offset + (i + 1) * ylen + j,
                                    vertex_offset + (i + 1) * ylen + j + 1,
                                    vertex_offset + i * ylen + j + 1,
                                ],
                                vertex_normal_indices=[],
                                texture_vertex_indices=[2, 1, 3],
                            )
                        )

    create_surface(True, 10.0)
    create_surface(False, 9.0, create_faces=False)
    create_surface(False, 0.0)

    def create_sides(top_layer_offset=0,bottom_layer_offset=xlen*ylen):
        # Top and bottom
        for i in range(0, ylen - 1):
            faces.append(
                Face(
                    vertex_indices=[
                        top_layer_offset+i,
                        top_layer_offset+i + 1,
                        bottom_layer_offset + i,
                    ],
                    vertex_normal_indices=[],
                    texture_vertex_indices=[0, 1, 2],
                )
            )
            faces.append(
                Face(
                    vertex_indices=[
                        bottom_layer_offset + i + 1,
                        bottom_layer_offset + i,
                        top_layer_offset+i + 1,
                    ],
                    vertex_normal_indices=[],
                    texture_vertex_indices=[2, 3, 0],
                )
            )
            # Top Side
            faces.append(
                Face(
                    vertex_indices=[
                        top_layer_offset+((xlen - 1) * ylen) + i + 1,
                        top_layer_offset+((xlen - 1) * ylen) + i,
                        (bottom_layer_offset) + (((xlen - 1) * ylen) + i),
                    ],
                    vertex_normal_indices=[],
                    texture_vertex_indices=[0, 1, 2],
                )
            )
            faces.append(
                Face(
                    vertex_indices=[
                        top_layer_offset+((xlen - 1) * ylen) + i + 1,
                        (bottom_layer_offset) + (((xlen - 1) * ylen) + i),
                        (bottom_layer_offset) + (((xlen - 1) * ylen) + i + 1),
                    ],
                    vertex_normal_indices=[],
                    texture_vertex_indices=[0, 1, 2],
                )
            )
        for i in range(0, xlen - 1):
            # Left Side
            faces.append(
                Face(
                    vertex_indices=[
                        top_layer_offset+(i + 1) * ylen,
                        top_layer_offset+i * ylen,
                        bottom_layer_offset + i * ylen,
                    ],
                    vertex_normal_indices=[],
                    texture_vertex_indices=[0, 1, 2],
                )
            )
            faces.append(
                Face(
                    vertex_indices=[
                        bottom_layer_offset + i * ylen,
                        bottom_layer_offset + (i + 1) * ylen,
                        top_layer_offset+(i + 1) * ylen,
                    ],
                    vertex_normal_indices=[],
                    texture_vertex_indices=[2, 3, 0],
                )
            )
            # Right Side
            faces.append(
                Face(
                    vertex_indices=[
                        top_layer_offset+(ylen - 1) + i * ylen,
                        top_layer_offset+(ylen - 1) + (i + 1) * ylen,
                        bottom_layer_offset + i * ylen + (ylen - 1),
                    ],
                    vertex_normal_indices=[],
                    texture_vertex_indices=[0, 1, 2],
                )
            )
            faces.append(
                Face(
                    vertex_indices=[
                        top_layer_offset+(ylen - 1) + (i + 1) * ylen,
                        bottom_layer_offset + (ylen - 1) + (i + 1) * ylen,
                        bottom_layer_offset + i * ylen + (ylen - 1),
                    ],
                    vertex_normal_indices=[],
                    texture_vertex_indices=[0, 1, 2],
                )
            )

    create_sides(0,xlen*ylen)
    create_sides(xlen*ylen,2*xlen*ylen)

    material = Material(colour=RGB(r=255.0, g=255.0, b=255.0))
    return Object(
        name="surface",
        material=material,
        verticies=verticies,
        faces=faces,
        vertex_normals=[],
        texture_verticies=texture,
    )


def write_object(obj: Object, filename: str):
    obj_name = obj.name.removesuffix(".obj")
    with open(filename, "w") as f:
        f.write("# Generated by uv script\n")
        f.write(f"o {obj_name}\n")
        f.write(f"mtllib {obj_name}.mtl\n")
        for v in obj.verticies:
            if isinstance(v, Point3):
                f.write(f"v {v.x} {v.y} {v.z}\n")
            elif isinstance(v, Point4):
                f.write(f"v {v.x} {v.y} {v.z} {v.w}\n")
            elif isinstance(v, Point6):
                f.write(f"v {v.x} {v.y} {v.z} {v.r} {v.g} {v.b}\n")
        for n in obj.vertex_normals:
            f.write(f"vn {n.x} {n.y} {n.z}\n")
        for t in obj.texture_verticies:
            if t.v is not None and t.w is not None:
                f.write(f"vt {t.u} {t.v} {t.w}\n")
            else:
                f.write(f"vt {t.u} {t.v or 0}\n")

        f.write(f"g {obj_name}\n")
        f.write(f"usemtl {obj_name}\n")
        print(len(obj.faces), "faces")
        for face in obj.faces:
            face.texture_vertex_indices = []
            if (
                len(face.vertex_normal_indices) == 0
                and len(face.texture_vertex_indices) == 0
            ):
                # If no normals or texture coordinates, just use vertex indices
                indices = " ".join([f"{o+1}" for o in face.vertex_indices])
                f.write(f"f {indices}\n")
                continue
            if len(face.texture_vertex_indices) == 0:
                indices = " ".join(
                    [
                        f"{o[0]+1}//{o[1]+1}"
                        for o in zip(face.vertex_indices, face.vertex_normal_indices)
                    ]
                )
                f.write(f"f {indices}\n")
                continue
            if len(face.vertex_normal_indices) == 0:
                indices = " ".join(
                    [
                        f"{o[0]+1}/{o[1]+1}"
                        for o in zip(face.vertex_indices, face.texture_vertex_indices)
                    ]
                )
                f.write(f"f {indices}\n")
                continue
            # Remaining valid condition is all are present
            indices = " ".join(
                [
                    f"{o[0]+1}/{o[1]+1}/{o[2]+1}"
                    for o in zip(
                        face.vertex_indices,
                        face.texture_vertex_indices,
                        face.vertex_normal_indices,
                    )
                ]
            )
            f.write(f"f {indices}\n")
            continue
    with open(filename.replace(".obj", ".mtl"), "w") as f:
        f.write(f"newmtl {obj_name}\n")
        colour = obj.material.colour
        f.write(f"  Kd {colour.r} {colour.g} {colour.b}")

def run():
    write_object(create_surface(), "surface.obj")

cProfile.run('run()', 'mesh-writer-1000.prof')
