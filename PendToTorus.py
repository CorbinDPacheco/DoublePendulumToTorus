from manim import *
from scipy.integrate import odeint 
import numpy as np

L1,L2,m1,m2 = 1,1,1,1 # lengths and masses of the pendula, respectively
    



class PendulumTorus(ThreeDScene):
    def construct(self):
        tmax, dt = 60, 0.05
        t = np.arange(0,tmax + dt, dt) 
        y0 = np.array([2.38,0,-1.82,0])
        y = odeint(self.deriv, y0, t, args = (L1,L2,m1,m2))
        theta1, theta2 = y[:,0], y[:,2]
        u, v = y0[0],y0[2]
        P = np.array([np.cos(u),np.sin(u),0])
        majorRadius = 3
        minorRadius = 1
        axes = ThreeDAxes()
        torus = Torus(major_radius = majorRadius, minor_radius = minorRadius) # u: big circle (azimuthal angle) , v: small circle
        point = Sphere(center = self.pointLocation(u,v,majorRadius,minorRadius,P), radius = 0.08).set_color(RED)
        self.set_camera_orientation(phi=50 * DEGREES, theta=45 * DEGREES)
        self.add(point,axes, torus)
        self.begin_ambient_camera_rotation(rate = 1)
        point.generate_target()
        for i in range(len(t)-1):
            P = np.array([np.cos(theta1[i+1]),np.sin(theta1[i+1]),0])
            pointF = Sphere(center = self.pointLocation(theta1[i+1],theta2[i+1],majorRadius,minorRadius,P), radius = 0.08).set_color(RED)
            self.play( Transform(point,pointF),run_time = 0.01)
        
    def pointLocation(self,u,v,majorRadius,minorRadius,P):
        return (majorRadius + minorRadius * np.cos(v)) * P + minorRadius * np.sin(v) * OUT
    
    def deriv(self, y, t, L1, L2, m1, m2):
        g = 9.8
        # Return the first derivatives of y = theta1, z1, theta2, z2.
        theta1, z1, theta2, z2 = y
        c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)
        theta1dot = z1
        z1dot = (m2*g*np.sin(theta2)*c - m2*s*(L1*z1**2*c + L2*z2**2) -
            (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
        theta2dot = z2
        z2dot = ((m1+m2)*(L1*z1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + 
            m2*L2*z2**2*s*c) / L2 / (m1 + m2*s**2)
        return theta1dot, z1dot, theta2dot, z2dot




class Pendula(Scene):
    def construct(self): 
        tmax, dt = 60, 0.05
        t = np.arange(0,tmax + dt, dt) 
        y0 = np.array([2.38,0,-1.82,0])
        y = odeint(self.deriv,y0,t,args = (L1,L2,m1,m2))
        theta1, theta2 = y[:,0], y[:,2]
        x1 = L1 * np.sin(theta1)
        y1 = -L1 * np.cos(theta1)
        x2 = x1 + L2 * np.sin(theta2)
        y2 = y1 - L2 * np.cos(theta2)
        
        line1 = Line(ORIGIN, [x1[0],y1[0],0])
        midJoint = Dot(point = [x1[0],y1[0],0], color = RED)
        line2 = Line([x1[0],y1[0],0],[x2[0],y2[0],0])
        endJoint = Dot(point = [x2[0],y2[0],0], color = RED)
        
        
        self.add(Dot(point = ORIGIN, color = RED),line1, line2, midJoint, endJoint)
        for i in range(len(t)-1):
            line1F = Line(ORIGIN, [x1[i+1],y1[i+1],0])
            line2F = Line([x1[i+1],y1[i+1],0],[x2[i+1],y2[i+1],0])
            midJointF = Dot(point = [x1[i+1],y1[i+1],0], color = RED)
            endJointF = Dot(point = [x2[i+1],y2[i+1],0], color = RED)
            
            self.play(
            Transform(line1,line1F),
            Transform(line2,line2F),
            Transform(midJoint,midJointF),
            Transform(endJoint,endJointF),
            run_time=0.01)
        
        
        
        
    def deriv(self, y, t, L1, L2, m1, m2):
        g = 9.8
        # Return the first derivatives of y = theta1, z1, theta2, z2.
        theta1, z1, theta2, z2 = y

        c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

        theta1dot = z1
        z1dot = (m2*g*np.sin(theta2)*c - m2*s*(L1*z1**2*c + L2*z2**2) -
                 (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
        theta2dot = z2
        z2dot = ((m1+m2)*(L1*z1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + 
                 m2*L2*z2**2*s*c) / L2 / (m1 + m2*s**2)
        return theta1dot, z1dot, theta2dot, z2dot

   

class Line3d(Scene):
    def construct(self):
        self.set_phi(45*DEGREES)
        self.set_theta(45*DEGREES)
        linex = Line([0,0,0],[1,0,0])
        liney = Line([0,0,0],[0,1,0])
        linez = Line([0,0,0],[0,0,1])
        self.play(GrowFromPoint(linex,ORIGIN))
        self.play(GrowFromPoint(liney,ORIGIN))
        self.play(GrowFromPoint(linez,ORIGIN))
        
