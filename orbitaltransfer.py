import pygame
import numpy
import math

import matplotlib.pyplot as plt

pygame.init()

WIDTH, HEIGHT = 800, 800
WIN = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Orbit simulation")

WHITE = (255, 255, 255)
YELLOW = (255, 255, 0)
GREEN = (0, 255, 0)
GREY = (192, 192, 192)
TEAL = (0, 128, 128)
CRIMSON = (220, 20, 60)
CYAN = (0, 255, 255)

top = 0
englist = []
ctlist = []
eclend = 0

def dist(point1, point2):
    xa, ya, za = point1
    xb, yb, zb = point2
    x = xb - xa
    y = yb - ya
    z = zb - za
    di = math.sqrt(x**2 + y**2 + z**2)
    return di

class Body:
    AU = 149.6e9
    ER =  6371000
    LEO = 2000000 + ER
    G = 6.67428e-11
    R_SUN = 696000000
    R_EARTH = 6371000
    SCALE = 250 / LEO
    TIMESTEP = 7
    SUN_LUM = 3.846e26
    EARTH_IRR = 240
    HOT_TIME = 0
    COLD_TIME = 0
    ECL_ANGLE = math.asin((R_SUN - R_EARTH)/AU)
    rho = 9.09e-17

    def __init__(self, x, y, z, radius, color, mass, real_rad):
        self.x = x
        self.y = y
        self.z = z
        self.radius = radius
        self.color = color
        self.mass = mass
        self.real_rad = real_rad

        self.orbit = []
        self.earth = False
        self.distance_to_earth = 0

        self.x_vel = 0
        self.y_vel = 0
        self.z_vel = 0
        if self.earth == False:
            self.apo = (0, 0, 0)
            self.per = (1e10, 1e10, 1e10)

            self.anode = (0, 0, 0)
            self.dnode = (0, 0 ,0)
            self.z_ant = 0
            self.per_ant = 0
            self.per_selected = True
            self.apo_selected = True
            self.anode_selected = False

            self.a = 0
            self.inc = 0
            self.ecc = 0
            self.arg = 0
            self.albedo = 0.1
            self.area = 36 * 280 + 14 * 100

            self.invol = 15 * math.pi * (140**2 - 125.8**2)
            self.airmass = self.invol * self.rho
            self.intemp = 0

            self.eng_cool = 2 * self.airmass * 1005 * (125-77) / 0.3
            self.eng_cool1 = self.eng_cool * 0.3
            self.cooler_pow = 0

    
    def toScale(self, point):
        x, y, z = point
        xn = x * self.SCALE + WIDTH/2
        yn = -y * self.SCALE + HEIGHT/2
        zn = -z * self.SCALE + HEIGHT/2
        return xn, yn, zn

    def elements(self, earth):
        apg = dist(self.apo, (earth.x, earth.y, earth.z))
        prg = dist(self.per, (earth.x, earth.y, earth.z))

        self.a = (apg + prg) / 2
        self.ecc = (apg - self.a) / self.a
        if self.per_selected == True and self.anode_selected == True:
            d1, d2, d3 = dist(self.anode, (self.per)), dist((self.per), (earth.x, earth.y, earth.z)), dist((earth.x, earth.y, earth.z), self.anode)
            cos_arg = (d2**2 + d3**2 - d1**2)/(2*d2*d3)
            self.arg = math.acos(cos_arg)



    def draw(self, win):
        global top
        x = self.x * self.SCALE + WIDTH/2
        y = -self.y * self.SCALE + HEIGHT/2
        z = -self.z * self.SCALE + HEIGHT/2
        if len(self.orbit) > 2:
            updated_points = []
            for var in self.orbit:
                x0, y0, z0 = self.toScale(var)    
                if top == 1:
                    updated_points.append((x0, z0))
                elif top == 0:
                    updated_points.append((x0, y0))

            pygame.draw.lines(win, self.color, False, updated_points, 4)
            if top == 0:
                pygame.draw.circle(win, self.color, (x, y), self.radius)
            elif top == 1:
                pygame.draw.circle(win, self.color, (x, z), self.radius)
        
    def drawElements(self, win, earth):
        self.elements(earth)
        global top
        ax, ay, az = self.toScale(self.apo)
        px, py, pz = self.toScale(self.per)


        x = earth.x*self.SCALE + WIDTH/2
        y = -earth.z*self.SCALE + HEIGHT/2
        z = -earth.z*self.SCALE + HEIGHT/2

        x_anode, y_anode, z_anode = self.toScale(self.anode)
        x_dnode, y_dnode, z_dnode = self.toScale(self.dnode)


        if top == 0:
            pygame.draw.line(win, TEAL, (x, y), (ax, ay) ,4)
            pygame.draw.line(win, CRIMSON, (x, y), (px, py), 4)
            pygame.draw.line(win, CYAN, (x_anode, y_anode), (x_dnode, y_dnode),4)
        elif top == 1:
            pygame.draw.line(win, TEAL, (x, z) , (ax, az), 4)
            pygame.draw.line(win, CRIMSON, (x, z), (px, pz), 4)
            pygame.draw.line(win, CYAN, (x_anode, z_anode), (x_dnode, z_dnode), 4)   
    
    def attraction(self, other):
        other_x, other_y, other_z = other.x, other.y, other.z
        distance_x = other_x - self.x
        distance_y = other_y - self.y
        distance_z = other_z - self.z
        distance = math.sqrt(distance_x**2 + distance_y**2 + distance_z**2)

        if  other.earth:
            self.distance_to_earth = distance

        force = self.G * self.mass * other.mass / distance**2
        theta = math.atan2(distance_y, distance_x)
        distance_a = math.sqrt(distance_x**2 + distance_y**2)
        psi = math.atan2(distance_z, distance_a)

        force_x = math.cos(psi) * math.cos(theta) * force
        force_y = math.cos(psi) * math.sin(theta) * force
        force_z = math.sin(psi) * force
        return force_x, force_y, force_z
    
    def drag(self):
        t_kelv = 273.1 -131.21 + self.distance_to_earth*0.00299
        pressure = 2.488 * (t_kelv/216.6)**(-11.388)
        rho = pressure/(0.2869 * t_kelv)

        self.airmass = self.invol * rho 
        vel = math.sqrt(self.x_vel**2 + self.y_vel**2 + self.z_vel**2)
        area = math.pi * self.real_rad**2
        drag = -0.5 * 0.5 * area * (vel**2) * rho

        psi = math.asin(self.z_vel/vel)
        theta = math.atan2(self.y_vel, self.x_vel)

        drag_plane = drag * math.cos(psi)
        drag_z = drag* math.sin(psi)
        drag_y = drag_plane * math.sin(theta)
        drag_x = drag_plane * math.cos(theta)

        return drag_x, drag_y, drag_z
    def heating(self):
        global eclend
        de = self.AU * self.R_EARTH/(self.R_SUN- self.R_EARTH)
        #print(de)
        if self.x <= de and self.x >=0 and abs(self.y) <= (de-self.x)* math.tan(self.ECL_ANGLE):
            self.COLD_TIME =self.COLD_TIME + self.TIMESTEP
            pow = self.area* (1-self.albedo) *  self.EARTH_IRR * (self.R_EARTH**2) / (dist((0, 0, 0), (self.x, self.y, self.z))**2)
            self.eng_cool = self.eng_cool +  2 * pow * self.TIMESTEP/0.3
            
            self.cooler_pow = self.eng_cool1/2149 + pow

            print(self.eng_cool)
            englist.append(self.eng_cool)
            ctlist.append(self.COLD_TIME)
            #print(ctlist)
        else:
            eclend = eclend + 1
            LUM = self.SUN_LUM * (1 - self.albedo) * self.area/ (4 * math.pi * (self.AU + self.x)**2)
            self.HOT_TIME = self.HOT_TIME + self.TIMESTEP
            #print("lum: " , LUM)
    def update_position(self, bodies):
        total_fx = total_fy = total_fz = 0
        for body in bodies:
            if self == body:
                continue
            fx, fy, fz = self.attraction(body)
            total_fx += fx
            total_fy += fy
            total_fz += fz
            if self.earth == False:
                drx, dry, drz = self.drag() #applying drag!!
                total_fx+= drx
                total_fy+= dry
                total_fz+= drz
        self.x_vel += total_fx / self.mass * self.TIMESTEP
        self.y_vel += total_fy / self.mass * self.TIMESTEP
        self.z_vel += total_fz / self.mass * self.TIMESTEP

        self.x += self.x_vel * self.TIMESTEP
        self.y += self.y_vel * self.TIMESTEP
        self.z += self.z_vel * self.TIMESTEP 
        self.orbit.append((self.x, self.y, self.z))

timestamp = 0
incr = 0
ecclist = []
dvlist = []
timestamps = []

def main():
    global top, incr, ecclist, dvlist, timestamp
    run = True
    clock = pygame.time.Clock()

    earth = Body(0, 0, 0, 50, GREEN, 5.7e24, Body.ER)
    earth.earth = True

    asteroid = Body(-1 * Body.LEO, 0,-10, 10 , GREY, 6000 , 500)
    asteroid.y_vel = 6753 + incr 
    bodies = [earth, asteroid]
    #top = int(input("0 for side, 1 for top "))
    print(Body.ECL_ANGLE * 180/ math.pi)
    while run:
        clock.tick(60)
        WIN.fill(WHITE)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False
        for body in bodies:
            if body.ecc >= 0.7:
                run = False
            body.z_ant = body.z
            body.per_ant = dist(body.per, (earth.x, earth.y, earth.z))
            body.update_position(bodies)
            timestamp = timestamp + 60 * body.TIMESTEP
            if body.earth == False:
                body.heating()
                if body.distance_to_earth > dist(body.apo, (earth.x, earth.y, earth.z)):
                    body.apo_selected = False
                    body.apo = (body.x, body.y, body.z)
                elif body.distance_to_earth < dist(body.apo, (earth.x, earth.y, earth.z)):
                    body.apo_selected= True
                if body.distance_to_earth < dist(body.per, (earth.x, earth.y, earth.z)):
                    body.per_selected = False
                    body.per = (body.x, body.y, body.z)
                elif body.distance_to_earth > dist(body.per, (earth.x, earth.y, earth.z)):
                    body.per_selected = True

                body.drawElements(WIN, earth)
                if body.z_ant * body.z <= 0 and body.z_ant < 0:
                    body.anode = (body.x, body.y, body.z)
                    body.anode_selected = True
                elif body.z_ant * body.z <= 0 and body.z_ant > 0:
                    body.dnode = (body.x, body.y, body.z)
                    body.anode_selected = False
                #if body.apo_selected == True and body.per_selected == True and dist(body.apo ,body.per) > body.a:
                    #ecclist.append(body.ecc)
                    #dvlist.append(incr)
                    #incr = incr+25
                    #timestamps.append(timestamp)
                    #print(ecclist)
                    #print(dvlist)
                    
            body.draw(WIN)

        pygame.display.update()
    e1 = 20312737834.023857
    e0 = 66551974.75528827

    print((e0)/asteroid.COLD_TIME)
    ifnal = 9421212.591562852
    intensity = 1000
    e_circ = ((2 * math.pi * 125.8 * 0.12276)**2) * intensity * asteroid.airmass * 1005 / (asteroid.cooler_pow * 0.3)
    en_tot = e1 +  2 * e_circ
    pygame.quit()


dvlist = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000, 1025, 1050, 1075, 1100, 1125, 1150, 1175, 1200, 1225, 1250, 1275, 1300, 1325, 1350, 1375, 1400, 1425, 1450, 1475, 1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725, 1750, 1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975, 2000, 2025]
ecclist = [0.0044566334755045825, 0.011254091161685954, 0.01858107538069753, 0.02601754164583963, 0.033508845627631414, 0.0410395993215805, 0.0486058036532956, 0.056202624542801276, 0.06382953680845038, 0.07148716624716277, 0.07917229319559202, 0.08688531685775196, 0.09462815041868569, 0.10239864347177, 0.11019648222287931, 0.11802217338062418, 0.12587541288658877, 0.1337557698848674, 0.1416647510536813, 0.14960340602284303, 0.1575679789269459, 0.16556040584683618, 0.17357968254517034, 0.18162775162381492, 0.18970190416932925, 0.1978051630464793, 0.20593466562379298, 0.214092331649774, 0.22227762897075104, 0.2304895215595021, 0.23873030100492287, 0.24699791517901493, 0.25529274396539586, 0.2636157145311236, 0.271965503204781, 0.2803429218787447, 0.28874834504507213, 0.2971815465226431, 0.30564142736175315, 0.31412915420602017, 0.32264482963893726, 0.33118782110858486, 0.3397631027026312, 0.3483615801466216, 0.3569872482862814, 
0.36564065960075104, 0.3743215957566877, 0.3830299806201872, 0.3917656871154426, 0.40052901381999667, 0.40932015712027964, 0.4181385680788858, 0.4269842940810532, 0.4358578565795065, 0.44475873205845146, 0.45368731686303526, 0.4626430843382386, 0.47162662632242836, 
0.48063758752656166, 0.4896760004782783, 0.49874220168712274, 0.5078357153845064, 0.5169567272084583, 0.5261052456838143, 0.5352811978254125, 0.5444846602976701, 0.5537158040070048, 0.5629743204581272, 0.5722604786447165, 0.5815739687078306, 0.5909151236672604, 0.6002836688172309, 0.60967977485946, 0.6191034041213307, 0.6285545120542387, 0.6380330926327904, 0.6475392374281811, 0.6570728642096819, 0.6666340217816309, 0.6762226521053574, 0.6858388110245311, 0.6954825176943648]


main()

plt.plot(ctlist, englist)

plt.title("Cooling Energy vs. Timestamp")
plt.ylabel("Cooling Energy (J)")
plt.xlabel("Timestamp (s)")


plt.show()