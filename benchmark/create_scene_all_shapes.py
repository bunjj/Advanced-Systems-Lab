from helper import *
from random import random as rand
import json
import sys

shapes = [ "box", "sphere", "cone" , "torus", "octahedron"]

if len(sys.argv) != 3:
    sys.exit("not the right number of arguments")

n = int(sys.argv[1])
path = sys.argv[2]

#create scene
scene = {}
scene["camera"] = camera()
scene["pointlight"] = pointlight()
objects = []
#load plane
o = object(kind = "plane")
o["params"] = plane()
objects.append(o)

for shape in shapes:
    for i in range(n):
        pos = [0,0,0]
        pos[0] = rand() * 16 - 8
        pos[1] = rand() * 6 - 3
        pos[2] = rand() * 5 + 12.5
        color = [rand(), rand(), rand()]
        rotation = [30*rand(),30*rand(),30*rand()]
        reflection = rand()

        o = object(kind = shape, position= pos, color= color, reflection=reflection, rotation=rotation)

        if shape == "box":
            extent = [rand(), rand(), rand()]
            o["params"] = box(extent)
        elif shape == "sphere":
            radius = rand() * 0.75
            o["params"] = sphere(radius)
        elif shape == "cone":
            o["params"] = cone([rand(), rand(), rand()])
        elif shape == "torus":
            r2 = rand()
            r1 = rand() + r2
            o["params"] = torus(r1, r2)
        elif shape == "octahedron":
            o["params"] = octahedron(rand())
        else:
            sys.exit("Wrong shape format")
        objects.append(o)


scene["objects"] =objects

with open(path, 'w') as outfile:
    json.dump(scene, outfile, indent= 4)
#print(json.dumps(scene, indent = 4))

