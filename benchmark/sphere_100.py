from helper import *
from random import random as rand
import json

scene = {}
scene["camera"] = camera()
scene["pointlight"] = pointlight()
objects = []
o = object(kind = "plane")
o["params"] = plane()
objects.append(o)

for i in range(200):
    pos = [0,0,0]
    pos[0] = rand() * 16 - 8
    pos[1] = rand() * 6 - 3
    pos[2] = rand() * 5 + 12.5
    color = [rand(), rand(), rand()]
    reflection = rand()
    radius = rand() * 0.75

    o = object(kind = "sphere", position= pos, color= color, reflection=reflection)
    o["params"] = sphere(radius)
    objects.append(o)


scene["objects"] =objects

with open('sphere_100.json', 'w') as outfile:
    json.dump(scene, outfile, indent= 4)
#print(json.dumps(scene, indent = 4))

