# This is the provided code for the ray tracing project.
#
# The most important part of this code is the command interpreter, which
# parses the scene description (.cli) files.
# Roger Lee

from __future__ import division
import traceback

debug_flag = False   # print debug information when this is True

    
def setup():
    size(320, 320) 
    noStroke()
    colorMode(RGB, 1.0)  # Processing color values will be in [0, 1]  (not 255)
    background(0, 0, 0)
    frameRate(30)

# make sure proper error messages get reported when handling key presses
def keyPressed():
    try:
        handleKeyPressed()
    except Exception:
        traceback.print_exc()

# read and interpret a scene description .cli file based on which key has been pressed
def handleKeyPressed():
    if key == '1':
        interpreter("01_one_sphere.cli")
    elif key == '2':
        interpreter("02_three_spheres.cli")
    elif key == '3':
        interpreter("03_shiny_sphere.cli")
    elif key == '4':
        interpreter("04_many_spheres.cli")
    elif key == '5':
        interpreter("05_one_triangle.cli")
    elif key == '6':
        interpreter("06_icosahedron_and_sphere.cli")
    elif key == '7':
        interpreter("07_colorful_lights.cli")
    elif key == '8':
        interpreter("08_reflective_sphere.cli")
    elif key == '9':
        interpreter("09_mirror_spheres.cli")
    elif key == '0':
        interpreter("10_reflections_in_reflections.cli")
    elif key == '-':
        interpreter("11_star.cli")

# You should add code for each command that calls routines that you write.
# Some of the commands will not be used until Part B of this project.
def interpreter(fname):
    global shapes
    global lights
    global eye
    global uvw
    global fov
    global bg_color
    reset_scene()  # you should initialize any data structures that you will use here
    
    fname = "data/" + fname
    # read in the lines of a file
    with open(fname) as f:
        lines = f.readlines()

    # parse the lines in the file in turn
    for line in lines:
        words = line.split()  # split up the line into individual tokens
        if len(words) == 0:   # skip empty lines
            continue
        if words[0] == 'sphere':
            x = float(words[2])
            y = float(words[3])
            z = float(words[4])
            radius = float(words[1])
            # call your sphere making routine here
            shapes.append(Sphere(radius,PVector(x,y,z),surface))
            # for example: create_sphere(x,y,z,radius)
        elif words[0] == 'fov':
            fov = float(words[1])
            pass
        elif words[0] == 'eye':
            eye = PVector(float(words[1]),float(words[2]),float(words[3]))
            pass
        elif words[0] == 'uvw':
            uvw.append(PVector(float(words[1]),float(words[2]),float(words[3])))
            uvw.append(PVector(float(words[4]),float(words[5]),float(words[6])))
            uvw.append(PVector(float(words[7]),float(words[8]),float(words[9])))
            pass
        elif words[0] == 'background':
            bg_color = (float(words[1]),float(words[2]),float(words[3]))
            pass
        elif words[0] == 'light':
            lights.append(Light(PVector(float(words[1]),float(words[2]),float(words[3])), PVector(float(words[4]),float(words[5]),float(words[6]))))
            pass
        elif words[0] == 'surface':
            surface = Material(PVector(float(words[1]),float(words[2]),float(words[3])), PVector(float(words[4]),float(words[5]),float(words[6])), PVector(float(words[7]),float(words[8]),float(words[9])), float(words[10]), (float(words[11])))
            pass
        elif words[0] == 'begin':
            aTriangle = Triangle([], surface) 
            pass
        elif words[0] == 'vertex':
            aTriangle.addVertex(PVector(float(words[1]),float(words[2]),float(words[3])))
            pass
        elif words[0] == 'end':
            shapes.append(aTriangle)
            pass
        elif words[0] == 'render':
            render_scene()    # render the scene (this is where most of the work happens)
        elif words[0] == '#':
            pass  # ignore lines that start with the comment symbol (pound-sign)
        else:
            print ("unknown command: " + word[0])

# render the ray tracing scene
def render_scene():
    global debug_flag

    for j in range(height):
        for i in range(width):
            debug_flag = False
            if i == 165 and j == 195:
                debug_flag = True
            
            d = (1/tan(radians(fov)/2.0))
            v = -((2*j)/height-1)
            u = (2*i)/width-1
            directionVector = PVector.normalize(PVector.add(PVector.add(PVector.mult(uvw[0],u),PVector.mult(uvw[1],v)), PVector.mult(uvw[2],-d)))
            ray = Ray(eye, directionVector)
            hit = rayIntersectScene(ray)
            # Maybe set a debug flag to true for ONE pixel.
            # Have routines (like ray/sphere intersection)print extra information if this flag is set.
            

            # create an eye ray for pixel (i,j) and cast it into the scene

            c = shadingEq(hit, ray)
            pix_color = color(c.x,c.y,c.z)
            set (i, j, pix_color)         # draw the pixel with the calculated color

# here you should reset any data structures that you will use for your scene (e.g. list of spheres)
def reset_scene():
    global shapes
    global lights
    global eye
    global uvw
    global fov
    global bg_color
    shapes = []
    lights = []
    uvw = []
    eye = 0
    fov = 0
    bg_color = 0
    pass

# prints mouse location clicks, for help debugging
def mousePressed():
    print ("You pressed the mouse at " + str(mouseX) + " " + str(mouseY))

# this function should remain empty for this assignment
def draw():
    pass
    
def rayIntersectScene(eye_ray):
    hit = None
    minT = 999999

    for shap in shapes:

        t = shap.intersect(eye_ray)

        if t != None and isinstance(shap, Triangle) and t < minT and t > 0:
            
            intersectionPoint = PVector(eye_ray.origin.x + t * eye_ray.direction.x, eye_ray.origin.y + t * eye_ray.direction.y, eye_ray.origin.z + t * eye_ray.direction.z)
            triple1 = PVector.dot(PVector.cross(PVector.sub(intersectionPoint, shap.vertices[0]), PVector.sub(shap.vertices[1], shap.vertices[0])), shap.normal)
            triple2 = PVector.dot(PVector.cross(PVector.sub(intersectionPoint, shap.vertices[1]), PVector.sub(shap.vertices[2], shap.vertices[1])), shap.normal)
            triple3 = PVector.dot(PVector.cross(PVector.sub(intersectionPoint, shap.vertices[2]), PVector.sub(shap.vertices[0], shap.vertices[2])), shap.normal)
            if ((triple1 > 0) == (triple2 > 0) == (triple3 > 0)):
                minT = t
                if PVector.dot(eye_ray.direction, shap.normal) < 0:

                    hit = Hit(shap, shap.normal, t, intersectionPoint)
                else:

                    hit = Hit(shap, PVector.mult(shap.normal,-1.0), t, intersectionPoint)


        elif (t > 0 and t < minT):
            minT = t
            intersectionPoint = PVector(eye_ray.origin.x + t * eye_ray.direction.x, eye_ray.origin.y + t * eye_ray.direction.y, eye_ray.origin.z + t * eye_ray.direction.z)
            hit = Hit(shap, PVector.normalize(PVector.sub(intersectionPoint,shap.position)), t, intersectionPoint)
    return hit
        
class Sphere:
    def __init__(self, radius, position, diffuse_color):
        self.radius = radius
        self.position = position
        self.diffuse_color = diffuse_color
        
    def intersect(self, ray):
        ux = ray.origin.x - self.position.x
        uy = ray.origin.y - self.position.y
        uz = ray.origin.z - self.position.z

        a = ray.direction.x**2+ray.direction.y**2 + ray.direction.z**2
        b = 2*(ray.direction.x*ux) +  2*(ray.direction.y*uy) +  2*(ray.direction.z*uz)
        c = ux**2 + uy**2 +uz**2 - self.radius**2

        if (b**2 - 4*a*c) < 0:
            return None
        sol1 = (-b - (b**2 - 4*a*c)**(1/2))/(2*a)
        sol2 = (-b + (b**2 - 4*a*c)**(1/2))/(2*a)

        if(sol1 > 0):
            if(sol2 > 0):
                return min(sol1,sol2)  
            else:
                return sol1
        if(sol2 > 0):
            return sol2 

class Triangle:
    def __init__(self, vertices, diffuse_color):
        self.vertices = vertices
        self.diffuse_color = diffuse_color
        
    def addVertex(self, vertex):
        self.vertices.append(vertex)
        if len(self.vertices) == 3:
            self.normal = PVector.cross(PVector.sub(self.vertices[1],self.vertices[0]), PVector.sub(self.vertices[2],self.vertices[1])).normalize()
    
    def intersect(self, ray):
        t = None
        denom = PVector.dot(self.normal, ray.direction)
        if denom != 0:
            t = PVector.dot(self.normal, PVector.sub(self.vertices[0], ray.origin)) / denom
        
        return t
        
class Light:
    def __init__(self, position, color):
        self.position = position
        self.col = color

class Hit:
    def __init__(self, shap, norm_vec, t, int_point):
        self.shap = shap
        self.norm_vec = norm_vec
        self.t = t
        self.int_point = int_point
        
class Material:
    def __init__(self, diffuse_rgb, ambient_rgb, specular_rgb, specular_power, k_refl):
        self.diffuse_rgb = diffuse_rgb
        self.ambient_rgb = ambient_rgb
        self.specular_rgb = specular_rgb
        self.specular_power = specular_power
        self.k_refl = k_refl
    def __str__(self): 
        return str( (self.diffuse_rgb,self.ambient_rgb,self.specular_rgb,self.specular_power,self.k_refl))
    
class Ray:
    def __init__(self, origin, direction):
        self.origin = origin
        self.direction = direction


def shadingEq(hit,ray, max_depth = 10):
    if (hit == None):
        return PVector(bg_color[0],bg_color[1],bg_color[2])
    #initialize 
    totalColor = PVector(0.0,0.0,0.0)
    if isinstance(hit.shap, Triangle):
        n_vector = hit.norm_vec
    else:
        n_vector = PVector.normalize(PVector.sub(hit.int_point, hit.shap.position))
        
    #reflection calc
    if max_depth > 0 and hit.shap.diffuse_color.k_refl > 0:
        reflRayOrigin = PVector.add(hit.int_point, PVector.mult(n_vector, 0.0001))
        scalar = PVector.dot(n_vector,PVector.mult(ray.direction, -1))*2
        reflRayDirection = PVector.add(PVector.mult(n_vector, scalar), ray.direction).normalize()
        reflRay = Ray(reflRayOrigin,reflRayDirection)
        reflHit = rayIntersectScene(reflRay)
        reflColor = PVector.mult(shadingEq(reflHit,reflRay,max_depth - 1), hit.shap.diffuse_color.k_refl)
        totalColor = PVector.add(totalColor, reflColor)
    for light in lights:
        color = PVector(0.0,0.0,0.0)
        shadowTerm = 1
        surface = hit.shap.diffuse_color

    
        #shadow calc
        shadowRayOrigin = PVector.add(hit.int_point, PVector.mult(n_vector, 0.0001))
        shadowRayDirection = PVector.sub(light.position, hit.int_point).normalize()
        shadowRay = Ray(shadowRayOrigin, shadowRayDirection)
        shadowHit = rayIntersectScene(shadowRay)
        
        if (shadowHit != None) and shadowHit.t < PVector.dist(light.position, hit.int_point):
            shadowTerm = 0
            pass

        #diffuse calc
        hitDiffuseMat = hit.shap.diffuse_color.diffuse_rgb 
        lightColor = light.col
        l_vector = PVector.normalize(PVector.sub(light.position, hit.int_point))
        
        diffuseCo = max(PVector.dot(n_vector,l_vector),0)
        color = PVector.add(PVector.mult(PVector.pairwise_mult(hitDiffuseMat,lightColor), diffuseCo),color)
        
        #spec calc
        h = PVector.sub(PVector.sub(light.position, hit.int_point), ray.direction).normalize()
        specularCo = max(PVector.dot(h, n_vector),0)**(surface.specular_power)
        specColor = PVector.mult(PVector.pairwise_mult(lightColor, surface.specular_rgb),specularCo)
        
        color = PVector.add(color,specColor)
        color = PVector.mult(color, shadowTerm)
        totalColor = PVector.add(totalColor, color)
        
    totalColor = PVector.add(totalColor, surface.ambient_rgb)
    return totalColor

class PVector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return str((self.x, self.y, self.z))

    @staticmethod
    def add(a, b):
        return PVector(
            a.x + b.x,
            a.y + b.y,
            a.z + b.z,
        )

    @staticmethod
    def sub(a, b):
        return PVector(
            a.x - b.x,
            a.y - b.y,
            a.z - b.z,
        )
        
    def mag(self):
        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def magSq(self):
        return self.x * self.x + self.y * self.y + self.z * self.z

    def copy(self):
        return PVector(self.x, self.y, self.z)

    def div(self, n):
        return PVector(
            a.x / n,
            a.y / n,
            a.z / n,
        )

    @staticmethod
    def dist(a, b):
        return PVector.sub(a, b).mag()
    
    @staticmethod
    def mult(a, n):
        return PVector(
            n * a.x,
            n * a.y,
            n * a.z,
        )

    @staticmethod
    def pairwise_mult(a, b):
        return PVector(
            a.x * b.x,
            a.y * b.y,
            a.z * b.z,
        )

    @staticmethod
    def dot(a, b):
        return a.x * b.x + a.y * b.y + a.z * b.z

    @staticmethod
    def cross(a, b):
        return PVector(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x,
        )

    @staticmethod
    def normalize(a):
        mag = sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
        return PVector(a.x / mag, a.y / mag, a.z / mag)

    def normalize(self):
        mag = sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
        self.x /= mag
        self.y /= mag
        self.z /= mag
        return self
