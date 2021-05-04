def vector(vec = [0,0,0]):
    v = {}
    v["x"] = vec[0]
    v["y"] = vec[1]
    v["z"] = vec[2]
    return v
def camera(fov= 30, position = [0,0,0], rotation=[0,0,0]):
    cam = {}
    cam["fov"] = fov
    cam["position"] = vector(position)
    cam["rotation"] = vector(rotation)
    return cam

def pointlight(position = [0,100,0], emission = [200,200,200]):
    light = {}
    light["position"] = vector(position)
    light["emission"] = vector(emission)
    return light

#returns parameters which all objects have in common
def object(kind, position =[0,0,0], rotation = [0,0,0], reflection = 0, shininess = 15, color =[1,1,1]):
    o = {}
    o["kind"] = kind
    o["position"] = vector(position)
    o["rotation"] = vector(rotation)
    o["reflection"] = reflection
    o["shininess"] = shininess
    o["color"] = vector(color)

    return o

#the following functions return the individual params for a kind of objects
def plane(normal =[0,1,0], displacement = -3):
    p = {}
    p["normal"] = vector(normal)
    p["displacement"] = displacement
    return p

def box(extents = [0.25, 0.5, 1]):
    return {"extents" : vector(extents)}

def sphere(radius = 3):
    return {"radius" : radius}

def cone(params = [1, 0.5, 1]):
    return params

def torus(r1 = 1, r2 = 0.5):
    return {"r1": r1, "r2": r2}

def octahedron( s = 1):
    return {"s": s}

