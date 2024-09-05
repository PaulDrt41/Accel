import pygame
import numpy
import math
import time

import matplotlib.pyplot as plt

pygame.init()

WIDTH, HEIGHT = 800, 800
WIN = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Orbit simulation")

WHITE = (255, 255, 255)
BLACK = (0, 0 ,0)
YELLOW = (255, 255, 0)
GREEN = (0, 255, 0)
GREY = (192, 192, 192)
TEAL = (0, 128, 128)
CRIMSON = (220, 20, 60)
CYAN = (0, 255, 255)
ORANGE = (255, 140, 0)


NOM = 196 # number of coils * 3
timestamp = 0
en_tot = 0

def distance(p1, p2):
    x1, y1 = p1
    x2, y2 = p2

    d= math.sqrt((x1 - x2)**2  + (y1 - y2)**2)
    return d
class dip:
    PERM = 1.25663706e-6
    def __init__(self, divol, difil):
        self.divol = divol
        self.difil = difil

        self.dimom =  (self.difil * self.divol)/self.PERM
    def f2d(self, other, dist):
        d1 = self.dimom
        d2 = other.dimom

        f = (self.PERM * 6 * d1 * d2) / (4 * math.pi *(dist**4))
        return f
class payload:
    global NOM
    RMIC = 14.2
    H = 0.06469285
    R = H * math.cos(2 * math.pi/NOM) / (1 - math.cos(2 * math.pi/NOM)) + 0.001
    SCALE = 350/((RMIC + R))
    TIMESTEP = 0.001
    def __init__(self, x, y, radius, mass, omega, color):
        self.x = x
        self.y = y
        self.omega = omega
        self.linvel = 0
        self.prevvel = 0
        self.acc = 0
        self.alpha = 0
        self.pseudoalpha = 0
        self.radius = radius
        self.mass = mass
        self.color = color

        self.trajectory = []
        
        self.divol = 0.4 * 0.4 * 0.4
        self.difiel = 1.47

        self.dimom = (self.difiel * self.divol)/(1.25663706e-6)


    def toScale(self, point):
        x, y = point
        xn = x*self.SCALE + HEIGHT/2
        yn = -y*self.SCALE + WIDTH/2

        return xn, yn
    def draw(self, win):
        x_sc, y_sc = self.toScale(((self.R +self.H)*math.cos(self.pseudoalpha), (self.R +self.H)*math.sin(self.pseudoalpha)))

        #if len(self.trajectory) > 2:
            #updated_points = []
            #for pt in self.trajectory:
                #x0, y0 = self.toScale(pt)
                #updated_points.append((x0, y0))
            #pygame.draw.lines(win, BLACK, False, updated_points)

        pygame.draw.circle(win, self.color, (x_sc, y_sc), self.radius)


    def update_position(self, magnet1, magnet2):
        global timestamp, en_tot

        d1 = distance((self.x, self.y), (magnet1.x, magnet1.y))
        d2 = distance((self.x, self.y), (magnet2.x, magnet2.y))
        #print("dist: " , d1, "dist2: ", d2 )
        #print("m1: ", distance((0, 0), (self.x, self.y)) )
        #print(self.alpha)
        if d1 != 0: beta = math.asin((distance((0, 0), (self.x, self.y))) * math.sin(self.alpha) / d1)
        elif d1 == 0: beta = math.pi/2

        if d2 != 0:beta2 = math.asin((distance((0, 0), (self.x, self.y))) * math.sin(magnet2.alpha-self.alpha) / d2)
        elif d2 == 0: beta2 = math.pi/2
        #print("beta " , beta, " beta2: ", beta2)

        if beta > 0:
            force_front= magnet1.magforce(self, d1 + magnet1.len/2)
            force_back = magnet2.magforce(self, d2 + magnet2.len/2)

            theta = beta - self.alpha
            psi = beta2 - (magnet2.alpha - self.alpha)
            
            force_x1 = force_front * math.sin(theta)
            force_x2 = force_back * math.sin(psi)

            force_y1 = force_front * math.cos(theta)
            force_y2 = force_back * math.cos(psi)

            force_c = self.mass * self.omega * self.omega * (self.R + self.H)

            force_y = force_c + force_y2 - force_y1 
            #print(force_c, force_y1, force_y2)
            force_x = 3 * (force_x1 + force_x2)
            self.prevvel = self.linvel
            self.omega = self.omega + (force_x * self.TIMESTEP)/(self.mass * (self.R + self.H))
            self.linvel = self.omega * (self.R + self.H)

            self.acc = force_x / self.mass

            #print("omega: ", self.omega, " ")
            self.alpha = self.alpha - self.omega * self.TIMESTEP
            self.pseudoalpha = (self.pseudoalpha - self.omega * self.TIMESTEP)%(2 * math.pi)
            self.x = (self.R  + self.H)*(math.cos(self.alpha))
            self.y = (self.R  + self.H)*math.sin(self.alpha)

            #en_mag1 = magnet1.amps * (magnet1.resistance**2) * self.TIMESTEP
            #en_mag2 = magnet2.amps * (magnet2.resistance**2) * self.TIMESTEP

            #en_rn = 3 * (en_mag1 + en_mag2) / 3600000000
            #en_tot = en_tot + en_rn
            timestamp = timestamp + self.TIMESTEP
            self.trajectory.append((self.x, self.y))
    def current(self, magnet):
        resistance = magnet.rpl * self.R * 2 * math.pi



class solenoid:
    global NOM
    PERM = 1.25663706e-6
    def __init__(self, amps, turns, len, radius, x, y):
        self.amps = amps
        self.turns = turns
        self.len = len
        self.radius = radius

        self.rpl = 0.12276 #mo/m
        self.wlen = self.radius * 2 * math.pi * self.turns
        self.wradius = self.len/(2 * self.turns)

        self.resistance = self.rpl * self.wlen

        self.x = x
        self.y = y
        self.alpha = 2 * math.pi / NOM

        self.n = self.turns/self.len
        self.constt = self.amps * self.n * self.PERM / 2
    def draw(self, win):
        global NOM
        alpha = math.acos(payload.R/(payload.R + payload.H))

        for i in range(0, NOM):
            x, y = payload.toScale(payload, (payload.R * math.cos(i * alpha), payload.R * math.sin(i * alpha)))
            pygame.draw.rect(win, ORANGE, (x-2.5, y-2.5, 5, 5))
    def magfield(self, x):
        a = (self.len/2 - x)/(math.sqrt((x - self.len/2)**2 + self.radius**2))
        b = (self.len/2 + x)/(math.sqrt((x + self.len/2)**2 + self.radius**2))

        mf_x = self.constt * (a + b)

        return mf_x

    def magforce(self, dipole, x):
        dm = dipole.dimom
        a = self.constt * dm * (self.radius**2)

        b = 1/(math.sqrt(((x - self.len/2)**2 + self.radius**2)**3))
        c = 1/(math.sqrt(((x + self.len/2)**2 + self.radius**2)**3))

        forc = a * (b + c)

        return forc
incr = 0
dvlist = []
tlist = []
ilist = []
acclist = []

def main():
    run = True
    clock = pygame.time.Clock()
    
    global NOM, incr, timestamp
    timestamp = 0
    current_i = 1000 + incr
    #mag1 = dip(0.4 * 0.4 * 0.4, 1.47)
    #mag2 = dip(0.4 * 0.4 * 0.4, 1.47)

    #F = (mag1.f2d(mag2, payload.H))
    #rads = 6000 * 2000 * 2000 / F
    #print (rads)
    box = payload((payload.R + payload.H), 0, 5, 6000, 0, CRIMSON)
    coil = solenoid(current_i, 2500, 0.5, 0.5, payload.R, 0)
    coil2 = solenoid(current_i, 2500, 0.5, 0.5, payload.R * math.cos(2 * math.pi / NOM), payload.R * math.sin(2 * math.pi / NOM))
    print(box.R)
    #print(coil.magforce(box, coil.len/2 + 0.1))

    #neom = 2 * math.pi/(math.acos(rads/(rads + box.H)))
    #print(neom)

    while run:
        clock.tick(60)
        WIN.fill(WHITE)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False
        box.update_position(coil, coil2)
        if box.alpha <= 0:
            print(timestamp, " " , box.omega, "linvel: ", box.linvel)
            tlist.append(timestamp)
            acclist.append(box.acc)
            dvlist.append(box.linvel)
            box.alpha = math.acos(box.R/(box.R + box.H))
            box.x = (box.R + box.H) * math.cos(box.alpha)
            box.y = (box.R + box.H) * math.sin(box.alpha)
            if box.linvel >= 2000:
                run = False
                #ilist.append(current_i)
                #incr = incr + 20
                #print(acclist)
                #print(ilist)
                #if current_i >= 1500:
                    #run = False
                #else:
                    #main()
        pygame.draw.circle(WIN, BLACK, (WIDTH/2, HEIGHT/2), ((box.R + box.RMIC) * box.SCALE), 3)
        pygame.draw.circle(WIN, GREY, (WIDTH/2, HEIGHT/2), ((box.R) * box.SCALE), 3)
        box.draw(WIN)
        coil.draw(WIN)
        pygame.display.update()
    pygame.quit()   
        
main()

acclist1 = [32.2481189839015, 32.89350878728641, 33.540137030499146, 34.18366032164061, 34.82829647858317, 35.47614343352399, 36.11799142243961, 36.76588943093503, 37.41149200250363, 38.05295480929604, 38.69770605941408, 39.34411242514489, 39.988802778511506, 40.63273235239119, 41.2808545280004, 41.92856295429003, 42.56870359222824, 43.216828387468034, 43.859697093339285, 44.5045600763556, 45.14986618627606, 45.79221747288718, 46.4420848302064, 47.08234359295479, 47.71997967619912, 48.3772660734024]
ilist1 = [1000, 1020, 1040, 1060, 1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320, 1340, 1360, 1380, 1400, 1420, 1440, 1460, 1480, 1500]

#plt.plot(tlist, dvlist)
#plt.title("Delta-V vs. Timestamp")
#plt.xlabel("Timestamp (s)")
#plt.ylabel("Delta-V (m/s)")


#other ideas:
#plt.plot(enlist, dvlist)
#plt.title("Delta-V vs. Energy")
#plt.ylabel("Delta-V (m/s)")
#plt.xlabel("Energy (MWh)")

plt.plot(tlist, acclist)
plt.title("Acceleration vs. Timestamp")
plt.xlabel("Timestamp (s)")
plt.ylabel("Acceleration (m/s^2)")

#plt.title("Delta-V vs. ")
#plt.xlabel("Delta-V")
#plt.ylabel("")

plt.show()