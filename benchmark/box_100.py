from helper import *
from random import random as rand
import json

#create scene
scene = {}
scene["camera"] = camera()
scene["pointlight"] = pointlight()
objects = []
#load plane
o = object(kind = "plane")
o["params"] = plane()
objects.append(o)

for i in range(200):
    pos = [0,0,0]
    pos[0] = rand() * 16 - 8
    pos[1] = rand() * 6 - 3
    pos[2] = rand() * 5 + 12.5
    color = [rand(), rand(), rand()]
    rotation = [30*rand(),30*rand(),30*rand()]
    reflection = rand()
    extent = [rand(), rand(), rand()]

    o = object(kind = "box", position= pos, color= color, reflection=reflection, rotation=rotation)
    o["params"] = box(extent)
    objects.append(o)


scene["objects"] =objects

with open('box_100.json', 'w') as outfile:
    json.dump(scene, outfile, indent= 4)
#print(json.dumps(scene, indent = 4))

