#!/usr/bin/env python
"""
Function to generate radial turbine stator. 
Author: Ingo Jahn 16/03/2015
V5 : - seperate the TOP BLK as two seperate block, which allows a smoother transitation for the stator inlet section.

"""



import os as os
import numpy as np
import scipy as sci
import shutil as sh
from getopt import getopt
import sys as sys 
#import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

from libprep3 import *
from e3_block import *


# Helper Functions:

##
##############################################################################
##############################################################################
class supersonic_mesh():
    #
    def __init__(self,R_out,R_throat,R_in,R1,R2,throat,alpha,N,throat_out,bl_ratio,corner_flag=0 ,d_chamfer=0.,throat_out_ratio=1., theta_0=0.,Z=0.):
        print "Initialising supersonic_mesh()"
        self.R_out = R_out
        self.R_throat = R_throat
        self.R_in = R_in
        self.R1 = R1
        self.R2 = R2
        self.throat = throat
        self.alpha = alpha  #alpha_trailing
        self.N = N
        self.d_theta = - np.pi * 2. / N
        self.throat_out = throat_out
        self.bl_ratio = bl_ratio       
        self.corner_flag = corner_flag
        self.d_chamfer = d_chamfer #[m] define the radius of the chamfer
        self.throat_out_ratio = throat_out_ratio
        self.theta_0= theta_0
        self.Z = Z
        print "HERE"
        self.calc_fixed_nodes1()
        self.calc_flex_nodes1() 
        self.Spiral_Length(self.XR1,self.YR2)
    #
    def calc_fixed_nodes1(self):
        # calculate main construction nodes        
        # Nodes on outlet, inlet and throat center
        self.X = Node(self.R_out * np.cos(self.theta_0),
                    self.R_out * np.sin(self.theta_0),self.Z,label="X") # Define the X node
        temp = np.arcsin(self.R_out * np.sin(np.pi-self.alpha) / self.R_throat) # why use np.pi - alpha?
        theta_1 = self.alpha-temp
        temp = np.arcsin(self.R_out * np.sin(np.pi-self.alpha) / self.R_in)
        theta_2 = self.alpha-temp
        self.Xi = Node(self.R_in * np.cos(self.theta_0+theta_2), self.R_in * np.sin(self.theta_0+theta_2),self.Z, label="Xi")
        self.Xt = Node(self.R_throat * np.cos(self.theta_0+theta_1), self.R_throat * np.sin(self.theta_0+theta_1),self.Z, label="Xt")
        # create nodes rotated by 1 blade passage
        self.Y = self.X.copy(); self.Y.label = "Y"; self.Y.rotate_about_zaxis(self.d_theta)
        self.Yi = self.Xi.copy(); self.Yi.label = "Yi"; self.Yi.rotate_about_zaxis(self.d_theta)
        self.Yt = self.Xt.copy(); self.Yt.label = "Yt"; self.Yt.rotate_about_zaxis(self.d_theta) 
        return 0
    #
    def calc_flex_nodes1(self):
        # these nodes should be adjusted and moved during optimisation process
        # nodes on circles defining throat
        self.XR1 = Node(self.Xt.x + (0.5*self.throat+self.R1) *np.sin(self.theta_0+self.alpha), 
                    self.Xt.y - (0.5*self.throat+self.R1)*np.cos(self.theta_0+self.alpha),
                    self.Z, label="XR1")
        self.XR2 = Node(self.Xt.x - (0.5*self.throat+self.R2) *np.sin(self.theta_0+self.alpha), 
                    self.Xt.y + (0.5*self.throat+self.R2)*np.cos(self.theta_0+self.alpha),
                    self.Z, label="XR2")
        # nodes at throat
        self.Xt1 = Node(self.Xt.x + 0.5*self.throat *np.sin(self.theta_0+self.alpha), 
                  self.Xt.y - 0.5*self.throat *np.cos(self.theta_0+self.alpha), self.Z,label="a")
        self.Xt2 = Node(self.Xt.x - 0.5*self.throat *np.sin(self.theta_0+self.alpha), 
                self.Xt.y + 0.5*self.throat *np.cos(self.theta_0+self.alpha),self.Z, label="a")
        # nodes defining throat at outlet
        t1 = self.throat_out * self.throat_out_ratio / (self.throat_out_ratio +  1.) 
        t2 = self.throat_out * 1./self.throat_out_ratio / (self.throat_out_ratio +  1.)
        dtheta1 = np.arctan((t1 / np.cos(self.alpha)) / self.R_out)# change to arctan, JQ 19/11/2016
        dtheta2 = np.arctan(-(t2 / np.cos(self.alpha)) / self.R_out)
        self.XO1 = Node(self.R_out * np.cos(self.theta_0-dtheta1),
                    self.R_out * np.sin(self.theta_0-dtheta1),self.Z,label="XO1") #?
        self.XO2 = Node(self.R_out * np.cos(self.theta_0-dtheta2), 
                    self.R_out * np.sin(self.theta_0-dtheta2),self.Z,label="XO2") #?
       
        # create nodes rotated by 1 blade passage
        self.YR2 = self.XR2.copy(); self.YR2.label = "YR2"; self.YR2.rotate_about_zaxis(self.d_theta)
        self.Yt2 = self.Xt2.copy(); self.Yt2.label = "Yt2"; self.Yt2.rotate_about_zaxis(self.d_theta)
        self.YO2 = self.XO2.copy(); self.YO2.label = "YO2"; self.YO2.rotate_about_zaxis(self.d_theta)

        # nodes defining the chamfer of the nozzle outlet
        if self.corner_flag == 1:
            self.Xm = Node( self.X.x + np.cos(self.theta_0 + self.alpha)*self.d_chamfer,
                            self.X.y + np.sin(self.theta_0 + self.alpha)*self.d_chamfer, self.Z, label = "Xm")
            self.Ym = self.Xm.copy(); self.Ym.label = "Ym"; self.Ym.rotate_about_zaxis(self.d_theta)

        #########################################
        #AB = Line(Node(0.,0.,self.Z),Node(1.,1,self.Z))

        #CD = Line(Node(0.,1.,self.Z),Node(1.,0.,self.Z))
        #crossection = crossection_point(AB,CD,1.e-8)
        #point = crossection.get_point()
        #print "_____crossection Point : ---" , point




        return 0

    def calc_tuning_nodes1(self,Xs_ratio,I0_ratio,I1_ratio,ITOP_ratio):
        self.Xs = Node((1.-Xs_ratio)*self.Xt.x + Xs_ratio*self.Xi.x, 
                        (1.-Xs_ratio)*self.Xt.y + Xs_ratio*self.Xi.y, 
                        self.Z, label="Xs")
        self.Ys = self.Xs.copy(); self.Ys.label = "Ys"; self.Ys.rotate_about_zaxis(self.d_theta)
        #
        temp = Arc(self.Xi,self.Yi,Node(0.,0.,self.Z))
        self.I0 = Node(temp.eval(I0_ratio),label="I0")
        self.I1 = Node(temp.eval(I1_ratio),label="I1")
        self.ITOP = Node(temp.eval(ITOP_ratio),label="ITOP")
        #
        self.Xt1_BL = Node((1.-self.bl_ratio*2)*self.Xt1.x + self.bl_ratio*2*self.Xt.x, 
                        (1.-self.bl_ratio*2)*self.Xt1.y + self.bl_ratio*2*self.Xt.y, 
                        self.Z, label="Xt1_BL")
        self.Yt2_BL = Node((1.-self.bl_ratio*2)*self.Yt2.x + self.bl_ratio*2*self.Yt.x, 
                        (1.-self.bl_ratio*2)*self.Yt2.y + self.bl_ratio*2*self.Yt.y, 
                        self.Z, label="Yt2_BL")    
        #
        temp = Arc(self.XO1,self.X,Node(0.,0.,self.Z))
        self.XO1_BL = Node(temp.eval(self.bl_ratio/2),label="XO1_BL")
        temp = Arc(self.YO2,self.Y,Node(0.,0.,self.Z))
        self.YO2_BL = Node(temp.eval(self.bl_ratio/2),label="YO2_BL")


        return 0
    #
    def Line_XO1_Xt1(self,beta1,alpha1="None",L1=0.1,L2=0.1,Node_D="None",Node_E="None"):
        A = self.Xt1.copy(); A.label = "a"
        G = self.XO1.copy(); G.label = "g"
        Length = ( (A.x - G.x)**2 + (A.y-G.y)**2 )**0.5
        B = Node(self.XR1.x - self.R1 * np.sin(self.theta_0+self.alpha+beta1), 
                  self.XR1.y + self.R1 *np.cos(self.theta_0+self.alpha+beta1), self.Z,label="b")
        #print "A,B,beta1" , A,B,beta1*360./(2*np.pi)
        C = Node(B.x - L1*Length * np.cos(self.theta_0+self.alpha+beta1),
                B.y - L1*Length * np.sin(self.theta_0+self.alpha+beta1),self.Z,label="c")
        if alpha1 == "None":
            alpha1 = self.alpha
        F = Node(G.x + L2*Length * np.cos(self.theta_0+ alpha1),
                G.y + L2*Length * np.sin(self.theta_0+ alpha1),self.Z,label="f")
        # create intermediary nodes if not defined
        if Node_D == "None" and Node_E == "None":
            temp = Bezier([B,C,F,G],"",0.,1.,1)
            D = Node(temp.eval(0.4),label="d")
            E = Node(temp.eval(0.6),label="e")
        else:
            if Node_D == "None":
                temp = Bezier([B,C,E,F,G],"",0.,1.,1)
                D = Node(temp.eval(0.4),label="d")
            if Node_E == "None":
                temp = Bezier([B,C,D,F,G],"",0.,1.,1)
                E = Node(temp.eval(0.6),label="e")
        #combine lines for output
        # Note lines are going backwards
        print "------->",A,B
        Line1 = Arc(B,A,self.XR1)
        Line0 = Bezier([G,F,E,D,C,B,A],"",0.,1.,1)
        #Line0 = Line(G,B)
        temp = Line0 #Polyline([Line0,Line1],"",0.,1.,1)


        if self.corner_flag == 1:    
            #print "Here", Line0,Line1
            length_ratio = self.d_chamfer/Length
            #print "---------length ratio:" , length_ratio
            self.XO3 = Node(temp.eval(1.8*length_ratio),label='XO3')  
            self.XO4 = Node(temp.eval(0.70*length_ratio),label='XO4')

            #starting to create the XO3_BL and XO4_BL
            temp_Line0 = move_and_stretch( self.XO1_BL , self.Xt1_BL , temp).Line()

            TMP = distance_point(temp_Line0,Line1,1.0e-3)
            XO3_temp = Node(temp.eval(1.80001*length_ratio),"")  # create a little increment point on the line to calculate the local tangents
            k1_temp = TMP.get_v_tangent(self.XO3,XO3_temp)
            line_1temp = TMP.p_k_line(self.XO3,k1_temp)
            XO3_BL_position_ratio = TMP.lines_intersection(line_1temp) # find the crosection point between line_1temp and temp_Line0
            self.XO3_BL = Node(temp_Line0.eval(XO3_BL_position_ratio),"")

            XO4_temp = Node(temp.eval(0.700001*length_ratio),"")
            k2_temp = TMP.get_v_tangent(self.XO4,XO4_temp)
            line_2temp = TMP.p_k_line(self.XO4,k2_temp)
            XO4_BL_position_ratio = TMP.lines_intersection(line_2temp) 
            self.XO4_BL = Node(temp_Line0.eval(XO4_BL_position_ratio),"")


            #define the line between block L0 and M0:
            # L0M0
            # The control point is the mirror point of XO3 through the line vertical to XO3_XO3_BL
            A = self.XO3
            B = self.XO3_BL
            D = self.Xm
            k1 = (B.y - A.y)/(B.x - A.x)
            k2 = -1 / k1
            a = k2
            b = -1
            c = -a * B.x + B.y

            # http://stackoverflow.com/questions/3306838/algorithm-for-reflecting-a-point-across-a-line
            d = (A.x + (A.y - c )* a ) / (1 + a**2)
            C_x = 2*d - A.x
            C_y = 2*d*a - A.y + 2*c

            C = Node(C_x,C_y, self.Z,'')

            self.L0M0 = Bezier([B,C,D],"",0.,1.,1)  

        return temp
    #

    def Line_YO2_Yt1(self,beta2,alpha2="None",L1=0.1,L2=0.1,Node_D="None",Node_E="None"):# actually Line_YO2_Yt2 not Yt1
        A = self.Xt2.copy(); A.label = "a"
        G = self.XO2.copy(); G.label = "g"
        Length = ( (A.x - G.x)**2 + (A.y-G.y)**2 )**0.5
        B = Node(self.XR2.x + self.R2 * np.sin(self.theta_0+self.alpha-beta2), 
                  self.XR2.y - self.R2 *np.cos(self.theta_0+self.alpha-beta2), self.Z,label="b")
        C = Node(B.x - L1*Length * np.cos(self.theta_0+self.alpha-beta2),
                B.y - L1*Length * np.sin(self.theta_0+self.alpha-beta2),self.Z,label="c")
        if alpha2 == "None":
            alpha2 = self.alpha
        F = Node(G.x + L2*Length * np.cos(self.theta_0+ alpha2),
                G.y + L2*Length * np.sin(self.theta_0+ alpha2),self.Z,label="f")
        # create intermediary nodes if not defined
        if Node_D == "None" and Node_E == "None":
            temp = Bezier([B,C,F,G],"",0.,1.,1)
            D = Node(temp.eval(0.4),label="d")
            E = Node(temp.eval(0.6),label="e")
        else:
            if Node_D == "None":
                temp = Bezier([B,C,E,F,G],"",0.,1.,1)
                D = Node(temp.eval(0.4),label="d")
            if Node_E == "None":
                temp = Bezier([B,C,D,F,G],"",0.,1.,1)
                E = Node(temp.eval(0.6),label="e")
        #combine lines for output
        Line1 = Arc(B,A,self.XR2)
        Line0 = Bezier([G,F,E,D,C,B],"",0.,1.,1)
        temp = Polyline([Line0,Line1],"",0.,1.,1)   
        temp.rotate_about_zaxis(self.d_theta)    # This is line YO2_Yt2

        if self.corner_flag == 1:
            """

            calculate the chamfer parameters.
            """
            self.r_chamfer = 0.11e-3 #[m] define the radius of the chamfer

            r,t,z = cart2polar(self.XO1.x,self.XO1.y,self.XO1.z)
            r = r + self.r_chamfer
            x,y,z = polar2cart(r,t,z)
            XO1_temp = Node(x,y,z,"")

            r,t,z = cart2polar(self.YO2.x,self.YO2.y,self.XO2.z)
            r = r + self.r_chamfer
            x,y,z = polar2cart(r,t,z)
            YO2_temp = Node(x,y,z,"")

            line_XO1temp_YO2temp = Arc(XO1_temp,YO2_temp,Node(0., 0., self.Z))

            # in order to make the two lines same evaluation direction:
            temp_l1 = Arc(A,B,self.XR2)
            temp_l0 = Bezier([B,C,D,E,F,G],"",0.,1.,1)
            temp_t = Polyline([temp_l1,temp_l0],"",0.,1.,1)
            temp_t.rotate_about_zaxis(self.d_theta)
            #print temp_t.eval(1.), self.YO2

            TEST=distance_point(line_XO1temp_YO2temp,temp_t,self.r_chamfer)
            m_r,s_r = TEST.find_point()
            #print "ratio",m_r,s_r

            self.YOr1      = Node(temp_t.eval(s_r),label="YOr1")
            self.O_chamfer = Node(line_XO1temp_YO2temp.eval(m_r),label="O_chamfer")

            r,t,z = cart2polar(self.O_chamfer.x,self.O_chamfer.y,self.O_chamfer.z)
            r = r-self.r_chamfer
            x,y,z = polar2cart(r,t,z)
            self.YOr0 = Node(x,y,z,label = "YOr0")
            self.line_chamfer = Arc(self.YOr0,self.YOr1,self.O_chamfer)
            print "SUCCESS!",self.line_chamfer.eval(1.),self.YOr1


            self.Line_YO2BL_Yt2BL = move_and_stretch (self.YO2_BL , self.Yt2_BL , temp).Line()
            TEST2 = distance_point(self.Line_YO2BL_Yt2BL,temp,1.0e-3)
            YOr1temp = Node(temp_t.eval(s_r+0.0001),"")
            tangent_v = TEST2.get_v_tangent(self.YOr1,YOr1temp)
            line_l = TEST2.p_k_line(self.YOr1,tangent_v)
            YOr1_BL_position_ratio = TEST2.lines_intersection(line_l)
            self.YOr1_BL = Node(self.Line_YO2BL_Yt2BL.eval(YOr1_BL_position_ratio),"")
            #print self.YOr1_BL

            delta_r = ((self.YOr1_BL.x - self.YOr1.x)**2+(self.YOr1_BL.y-self.YOr1.y)**2)**0.5
            r,t,z = cart2polar(self.YOr0.x,self.YOr0.y,self.YOr0.z)
            r = r - delta_r
            x,y,z = polar2cart(r,t,z)
            self.YOr0_BL = Node(x,y,z,label = "YOr0_BL")
            self.Line_YOr0BL_YOr1_BL = Arc(self.YOr0_BL,self.YOr1_BL,self.O_chamfer)
            #print self.Line_YOr0BL_YOr1_BL.eval(0.),self.YOr0_BL

            r,t,z = cart2polar(self.XO1.x,self.XO1.y,self.XO1.z)
            r = r - delta_r
            #t = t + 0.4/360*2*np.pi
            x,y,z = polar2cart(r,t,z)
            self.XO1_BL_2 = Node(x,y,z,label = "XO1_BL_2")

            ##self.YOr1 = Node(temp.eval(self.ratio_chamfer_N),label="YOr1")
            F_y = F.copy();  F_y.rotate_about_zaxis(self.d_theta)
            E_y = E.copy();  E_y.rotate_about_zaxis(self.d_theta)
            D_y = D.copy();  D_y.rotate_about_zaxis(self.d_theta)
            C_y = C.copy();  C_y.rotate_about_zaxis(self.d_theta)
            B_y = B.copy();  B_y.rotate_about_zaxis(self.d_theta)
            A_y = A.copy();  A_y.rotate_about_zaxis(self.d_theta)
            temp_Line0= Bezier([self.YOr1,E_y,D_y,C_y,B_y],"",0.,1.,1)
            temp_Line1= Arc(B_y,A_y,self.YR2) 
    #        print "---------some points", self.YOr1,F_y,E_y,D_y,C_y,B_y,A_y
            self.Line_Yt2_YOr1 = Polyline([temp_Line0,temp_Line1],"",0.,1.,1)

            #temp_Line2 = move_and_stretch( self.YO2_BL , self.Yt2_BL , temp)
            #length_stretch = ((self.Yt2_BL.x - self.YO2_BL.x)**2 + (self.Yt2_BL.y - self.YO2_BL.y)**2)**0.5
            #length_ratio = self.d_chamfer*1.2/length_stretch
            #self.YOr1_BL = Node(temp_Line2.Line().eval(length_ratio), label = "YOr1_BL")

        return temp
    #
    #def Line_YOr0_YOr1(self):  
    #     r = self.d_chamfer / 6
    #     A = self.YOr0 # defined in calc_flex_nodes1
    #     B = self.YOr1 # defined in Line_YO2_Yt1
    #     C = Node( 0.5*(A.x + B.x), 0.5*(A.y + B.y), self.Z)
    #     q = ( (A.x - B.x)**2 + (A.y-B.y)**2 )**0.5        
         # chamfer center calculation
         # reference:  http://mathforum.org/library/drmath/view/53027.html
    #     O_chamfer_x     = C.x + ((r**2. - (q/2.)**2.)**0.5)*( A.y - B.y)/q
    #     O_chamfer_y     = C.y + ((r**2. - (q/2.)**2.)**0.5)*(B.x - A.x)/q
         #self.O_chamfer = Node( O_chamfer_x, O_chamfer_y, self.Z, label = "O_chamfer") 
         #Line = Arc(A,B,self.O_chamfer)
    #     return Line
    #
    def Path_Xt1_Yt2_1(self,t):   # not coming to here
        #print "Here2", t
        #PyFunctionPath(t) to evaluate line
        start = self.XR1.copy() # center of circle at start
        end = self.YR2.copy() # center of circle at end
        #
        L1 = ( (start.x-end.x)**2 + (start.y-end.y)**2 )**0.5
        #
        R_start = ( (self.Xt1.y - start.y)**2 + (self.Xt1.x - start.x)**2 )**0.5
        R_end   = ( (self.Yt2.y - end.y  )**2 + (self.Yt2.x - end.x  )**2 )**0.5
        #
        theta_start = np.arctan2(self.Xt1.y - start.y ,self.Xt1.x - start.x)
        theta_end   = np.arctan2(self.Yt2.y - end.y   ,self.Yt2.x - end.x  )
        theta_m     = np.arctan2(start.y-end.y         ,start.y-end.y) - np.pi/2. + np.arctan2(R_start-R_end,L1) #????
        #print 'different',np.arctan2(start.y-end.y         ,start.x-end.x), np.arctan2(start.y-end.y         ,start.y-end.y)
        # 
        L0 = R_start * (theta_start-theta_m)
        L2 = R_end * (theta_m - theta_end)
        LL = L0+L1+L2
        if t < L0/LL:
            tt = t * LL/L0
            C_x = start.x 
            C_y = start.y 
            R = R_start
            theta = (1.-tt)*theta_start + tt*theta_m            
        elif t > (L0+L1)/LL:
            tt = (t - (L0+L1)/LL) * LL/L2
            C_x = end.x 
            C_y = end.y 
            R = R_end
            theta = (1.-tt)*theta_m + tt*theta_end   
        else:
            tt = (t - L0/LL) * LL/L1
            C_x = (1.-tt)*start.x + tt*end.x
            C_y = (1.-tt)*start.y + tt*end.y 
            theta = theta_m #- np.pi/2.   
            R = (1.-tt)*R_start + tt*R_end        
        # 
        x = C_x + R * np.cos(theta)
        y = C_y + R * np.sin(theta)
        z = 0.0
        return (x,y,z)
    #
    def Path_Xt1_Yt2_2(self,t): # not coming here either
        #print "Here2", t
        #PyFunctionPath(t) to evaluate line
        start = self.XR1.copy() # center of circle at start
        end = self.YR2.copy() # center of circle at end
        #
        L1 = ( (start.x-end.x)**2 + (start.y-end.y)**2 )**0.5
        #
        R_start = ( (self.Xt1.y - start.y)**2 + (self.Xt1.x - start.x)**2 )**0.5
        R_end   = ( (self.Yt2.y - end.y  )**2 + (self.Yt2.x - end.x  )**2 )**0.5
        #
        theta_start = np.arctan2(self.Xt1.y - start.y ,self.Xt1.x - start.x)
        theta_end   = np.arctan2(self.Yt2.y - end.y   ,self.Yt2.x - end.x  )
        # Issues arise if theta_start < theta_end
        if theta_start < theta_end:
            theta_start = theta_start + 2.*np.pi
        theta_m     = np.arctan2(start.y-end.y         ,start.x-end.x) - np.pi/2.
        L0 = R_start * (theta_start-theta_m)
        L2 = R_end * (theta_m - theta_end)
        LL = L0+L1+L2
        if t < L0/LL:
            tt = t * LL/L0
            C_x = start.x 
            C_y = start.y 
            #R = R_start
            # R = (1.-tt)*R_start + tt*R_end  # linear variation
            # R = R_start + tt**2 * (R_end-R_start)   # rapid variation at end 
            # R = R_end + (1.-tt)**2 * (R_start-R_end)   # rapid variation at start
            R = R_start + (R_end-R_start) * (1. - np.cos(tt*np.pi)) *0.5
            theta = (1.-tt)*theta_start + tt*theta_m            
        elif t > (L0+L1)/LL:
            tt = (t - (L0+L1)/LL) * LL/L2
            C_x = end.x 
            C_y = end.y 
            R = R_end
            theta = (1.-tt)*theta_m + tt*theta_end   
        else:
            tt = (t - L0/LL) * LL/L1
            C_x = (1.-tt)*start.x + tt*end.x
            C_y = (1.-tt)*start.y + tt*end.y 
            theta = theta_m #- np.pi/2.  
            R = R_end 
        # 
        x = C_x + R * np.cos(theta)
        y = C_y + R * np.sin(theta)
        z = 0.0
        #print x,y,z
        return (x,y,z)
    ##
    def Line_XR1t_YR2t(self,X1s_ratio,TOP_ratio):
        
        #TOP_ratio = 0.25
        print "---------------------"
        print "TEST LINE_XR1T_YR2T"
        print "---------------------"

        self.nose_path = PyFunctionPath(self.Path_Xt1_Yt2_3)
        self.nose_path_bl= offset(self.Xt1_BL , self.Yt2_BL , self.nose_path).Line()
        # The following code is used to calculate XR1s, XR1s_BL0 and XR1s_BL
        self.XR1s     =  Node(self.nose_path.eval(X1s_ratio),label = "XR1s")
        self.TOP      =  Node(self.nose_path.eval(TOP_ratio),label = "TOP")
        self.XR1s_BL0 =  Node(self.nose_path_bl.eval(X1s_ratio),label = "XR1s_BL0")
        self.TOP_BL   =  Node(self.nose_path_bl.eval(TOP_ratio),label = "TOP_BL")
        self.L1_BL_S   = Polyline([self.nose_path],"",0.,X1s_ratio)
        self.TOPL_BL_W = Polyline([self.nose_path],"",TOP_ratio,X1s_ratio)

        #        d       nd
        #     o-----o----------o
        #     A     B          C
        #   
        #     A (Ax,Ay)   B(Bx,By)  C(Cx,Cy)
        #
        #   calculate the point XR1s_BL 

        A = self.XR1s
        B = self.XR1s_BL0
        
        n = 3.0   # set how far the point c should located.
        d = np.sqrt( (B.y - A.y)**2 + (B.x - A.x)**2 )
        tol = 1.e-8

        k = (B.y - A.y)/(B.x - A.x)

        cx1 = B.x - n* d /(np.sqrt(k**2+1))
        cx2 = B.x + n* d /(np.sqrt(k**2+1))

        cy1 = k*(cx1 - B.x) + B.y
        cy2 = k*(cx2 - B.x) + B.y

        d1 = np.sqrt( (cy1-A.y)**2 + (cx1 - A.x)**2 )
        d2 = np.sqrt( (cy2-A.y)**2 + (cx2 - A.x)**2 )
        if abs(d1 - (n+1)*d) < tol:
            cy = cy1
            cx = cx1
        elif abs(d2 - (n+1)*d) < tol:
            cy = cy2
            cx = cx2
        else:
            print "ERROR happens in calculating node XR1s_BL"
        self.XR1s_BL = Node(cx,cy, label = "XR1s_BL")


        # Then calculation line Xt1BL_XR1sBL
        self.Xt1BL_XR1sBL = move_and_stretch(self.Xt1_BL,self.XR1s_BL,self.L1_BL_S).Line() 
        self.TOPBL_XR1SBL = move_and_stretch(self.TOP_BL,self.XR1s_BL,self.TOPL_BL_W).Line() 

        return 
    ##
    def spiral_Line(self,r_start,r_end,t_start,t_end,t):
        # function to calculate position along spiral
        #print "r_start",r_start
        R = (1.-t)*r_start + t*r_end
        T = (1.-t)*t_start + t*t_end
        X = R*np.cos(T); Y = R*np.sin(T)
        return X,Y      
    ##
    def spiral_angle(self,r_start,r_end,t_start,t_end,t,dt=0.001):
        # function to get angle of sprial line
        t_p = t+dt
        if t_p > 1.: 
            t_p =1.
        t_m = t-dt
        if t_m < 0.:
            t_m = 0.
        #DT = t_p - t_m
        X0,Y0 = self.spiral_Line(r_start,r_end,t_start,t_end,t_p)
        X1,Y1 = self.spiral_Line(r_start,r_end,t_start,t_end,t_m)
        return np.arctan2(Y0-Y1,X0-X1)
    ##
    def Spiral_Length(self,start,end):
        #
        r_start = np.sqrt(start.x**2 + start.y**2)
        t_start = np.arctan2(start.y,start.x)
        r_end = np.sqrt(end.x**2 + end.y**2)
        t_end = np.arctan2(end.y,end.x)
        # 
        L1 = 0.
        t_list = np.arange(100)/99.
        for i in range(99):
            x0,y0 = self.spiral_Line(r_start,r_end,t_start,t_end,t_list[i])
            x1,y1 = self.spiral_Line(r_start,r_end,t_start,t_end,t_list[i+1])          
            L1 = L1 + np.sqrt( (x0-x1)**2 + (y0-y1)**2) 
        self.L1 = L1
        return   
    ##
    def Path_Xt1_Yt2_3(self,t):
        #print "Here3", t
        #PyFunctionPath(t) to evaluate line
        start = self.XR1.copy() # center of circle at start
        end = self.YR2.copy() # center of circle at end
        #
        r_start = np.sqrt(start.x**2 + start.y**2)
        t_start = np.arctan2(start.y,start.x)
        r_end = np.sqrt(end.x**2 + end.y**2)
        t_end = np.arctan2(end.y,end.x)
        #
        R_start = ( (self.Xt1.y - start.y)**2 + (self.Xt1.x - start.x)**2 )**0.5
        R_end   = ( (self.Yt2.y - end.y  )**2 + (self.Yt2.x - end.x  )**2 )**0.5
        #
        theta_start = np.arctan2(self.Xt1.y - start.y ,self.Xt1.x - start.x)
        theta_end   = np.arctan2(self.Yt2.y - end.y   ,self.Yt2.x - end.x  )
        #print "Angles1",theta_start,theta_end
        # Issues arise if theta_start < theta_end
        #if theta_start < theta_end:
        #    theta_start = theta_start + 2.*np.pi
        # theta_m     = np.arctan2(start.y-end.y         ,start.x-end.x) - np.pi/2.
        theta_mstart = self.spiral_angle(r_start,r_end,t_start,t_end,0.)+np.pi/2.
        theta_mend = self.spiral_angle(r_start,r_end,t_start,t_end,1.)+np.pi/2.

        """
        print "theta_mend", theta_mend
        print "theta_end", theta_end
        print "Diff", abs(theta_mend - theta_end)
        print R_start, R_end
        print "Done"
        """

        #print "Points:", r_start,t_start,r_end,t_end
        #print "Angles:", theta_start, theta_mstart, theta_mend, theta_end
        #print "delta_angle", theta_start-theta_mstart
        if theta_start < theta_mstart:
            theta_start += 2.*np.pi
        #print "Angles after:", theta_start, theta_mstart, theta_mend, theta_end
        L0 = R_start * (theta_start-theta_mstart)
        L1 = self.L1
        if theta_mend >= theta_end:
            L2 = R_end * (theta_mend - theta_end)
            LL = L0+L1+L2
            #print L0, L1, L2, LL
            if t < L0/LL:
                tt = t * LL/L0
                C_x = start.x 
                C_y = start.y 
                #R = R_start
                # R = (1.-tt)*R_start + tt*R_end  # linear variation
                # R = R_start + tt**2 * (R_end-R_start)   # rapid variation at end 
                # R = R_end + (1.-tt)**2 * (R_start-R_end)   # rapid variation at start
                R = R_start + (R_end-R_start) * (1. - np.cos(tt*np.pi)) *0.5
                theta = (1.-tt)*theta_start + tt*theta_mstart            
            elif t > (L0+L1)/LL:
                tt = (t - (L0+L1)/LL) * LL/L2
                C_x = end.x 
                C_y = end.y 
                R = R_end
                theta = (1.-tt)*theta_mend + tt*theta_end   
            else:
                tt = (t - L0/LL) * LL/L1
                C_x,C_y = self.spiral_Line(r_start,r_end,t_start,t_end,tt)
                theta = self.spiral_angle(r_start,r_end,t_start,t_end,tt)+np.pi/2.
                R = R_end 
        else:
            L2 = 0. # no part at end.
            LL = L0+L1+L2
            #print L0, L1, L2, LL
            if t < L0/LL:
                tt = t * LL/L0
                C_x = start.x 
                C_y = start.y 
                #R = R_start
                # R = (1.-tt)*R_start + tt*R_end  # linear variation
                # R = R_start + tt**2 * (R_end-R_start)   # rapid variation at end 
                # R = R_end + (1.-tt)**2 * (R_start-R_end)   # rapid variation at start
                R = R_start + (R_end-R_start) * (1. - np.cos(tt*np.pi)) *0.5
                theta = (1.-tt)*theta_start + tt*theta_mstart  
            else:
                tt = (t - L0/LL) * LL/L1
                theta_correction = tt * (theta_mend - theta_end)
                C_x,C_y = self.spiral_Line(r_start,r_end,t_start,t_end,tt)
                theta = self.spiral_angle(r_start,r_end,t_start,t_end,tt)+np.pi/2. - theta_correction
                R = R_end 
        # 
        x = C_x + R * np.cos(theta)
        y = C_y + R * np.sin(theta)
        z = self.Z
        # print "Inseide function:",t,tt, x,y,z
        return (x,y,z)

    def Line_Xm_XO3(self):
        XO3 = self.XO3
        U4  = self.Xm
        B8  = self.XO3_BL
        Line0 = Line(XO3,B8)
        Line1 = Line(B8,U4)
        temp = Polyline([Line0,Line1],"",0.,1.)
        #print "I'm here"
        #temp = Arc3(XO3,B8,U4)  # the last one is "1" rather than "1.", a string        
        return temp

    def Initialize_Bezier_Curves(self):


        # Yt2bl_Yt2:
        A = self.Yt2_BL
        B = self.Yt2
        AB = Line(A,B)
        C = Node(AB.eval(0.25),'')
        D = Node(AB.eval(0.5),'')
        E = Node(AB.eval(0.75),'')
        self.Yt2bl_Yt2 = Bezier([A,C,D,E,B],"",0.,1.,1)


        # Yt_Yt2bl
        A = self.Yt
        B = self.Yt2_BL
        AB = Line(A,B)
        C = Node(AB.eval(0.25),'')
        D = Node(AB.eval(0.5),'')
        E = Node(AB.eval(0.75),'')   
        self.Yt_Yt2bl =   Bezier([A,C,D,E,B],"",0.,1.,1)   


        # U0M1 
        # The control point is the mirror point of YOr1_BL through the line vertical to YOr1_YOr1_BL
        A = self.YOr1
        B = self.YOr1_BL
        D = self.Ym
        k1 = (B.y - A.y)/(B.x - A.x)
        k2 = -1 / k1
        a = k2
        b = -1
        c = -a * B.x + B.y

        # http://stackoverflow.com/questions/3306838/algorithm-for-reflecting-a-point-across-a-line
        d = (A.x + (A.y - c )* a ) / (1 + a**2)
        C_x = 2*d - A.x
        C_y = 2*d*a - A.y + 2*c

        C = Node(C_x,C_y, self.Z,'')

        self.U0M1 = Bezier([D,C,B],"",0.,1.,1)

        return






    
##
##############################################################################
##############################################################################
def get_ratio(path,point):  
    """
     This function is used to get the ratio for a point on a specific line
         path  - which path has the point
         point - the point will be get with the specific evaluation ratio
     Jianhui Qi
    """     
    high  = 1.
    low   = 0.

    for i in range (1000):
        mid  = 0.5 * (high + low)
        if (point.y - path.eval(mid).y)*(path.eval(high).y - path.eval(mid).y) < 0:
            high = mid
        elif (point.y - path.eval(mid).y)*(path.eval(high).y - path.eval(mid).y) > 0 :
            low = mid
        else:
            ratio = mid
            #print "\n Point", point, "ratio", ratio
            #print path.eval(ratio)
            break
        i += 1
            #print "i=",i
        
    if np.abs(path.eval(ratio).x - point.x) >= 1.e-7 :
        high  = 1.
        low   = 0.

        for j in range (1000):
             mid  = 0.5 * (high + low)
             if (point.x - path.eval(mid).x)*(path.eval(high).x - path.eval(mid).x) < 0:
                 high = mid
             elif (point.x - path.eval(mid).x)*(path.eval(high).x - path.eval(mid).x) > 0 :
                 low = mid
             else:
                 ratio = mid
                 #print "\n Point", point, "ratio", ratio
                 #print path.eval(ratio)
                 break
             j += 1
             #print "j=",j
        if np.abs(path.eval(ratio).y - point.y) >= 1.e-7:
             return "The evaluation of ratio for a specific point on a specific line is failed",1/0
    return ratio
  

##
##############################################################################
##############################################################################
class crossection_point():
    """
    calss that return the crossection point of two given lines.
        master_Path - path 1, the return point is on master path
        slave_Path  - path 2
        d           - the distance between return point and the slave path.
                    - if d > 1e-8, then return the point on master path which
                      has a distance of d.
    Author: Jianhui Qi 23/11/2016
    """
    def __init__(self,master_Path,slave_Path,d):

        self.AB = master_Path
        self.CD = slave_Path
        self.r_t = 1.e-10      # ratio tolerance, if the increment of path evaluation
                              # less than r_t, then the sign of delta_i change to 
                              # negative
        self.d_t = d
        self.init_S_i = 1.         # sign of delta_i
        self.init_S_j = 1.         # sign of delta_j
        self.i   = 0.         # initial eval ratio of master path
        self.j   = 0.         # initial eval ratio of slave path
        self.init_delta_i = 0.01   # initial increment of master path evaluation ratio
        self.init_delta_j = 0.01   # initial increment of salve path evaluation ratio
        self.A0 = Node(self.AB.eval(self.i),"")
        self.C0 = Node(self.CD.eval(self.i),"")
        self.max_iter = 400000     # define the maximum iteration number 
        self.S_j = self.init_S_j 
        self.S_i = self.init_S_i
        self.delta_i = self.init_delta_i
        self.delta_j = self.init_delta_j  
    #
    def find_A0(self):
        print "________________I'm in A0"
        #self.S_i = self.init_S_i
        #
        #self.A0 = self.AB.eval(self.i)
        self.A1 = Node(self.AB.eval(self.i + self.delta_i),"")
        A0C0 =  ((self.A0.x - self.C0.x)**2 + (self.A0.y - self.C0.y)**2 )**0.5
        A1C0 =  ((self.A1.x - self.C0.x)**2 + (self.A1.y - self.C0.y)**2 )**0.5

        while A1C0 - A0C0 > 1e-10:
            print "________________I'm in A02"
            print A1C0-A0C0, self.delta_i
            self.delta_i = 0.5 * self.delta_i
            self.A1 = Node(self.AB.eval(self.i + self.delta_i))
            if np.abs(self.delta_i) < self.r_t:
                self.delta_i = -self.init_delta_i + self.delta_i
            A1C0 =  ((self.A1.x - self.C0.x)**2 + (self.A1.y - self.C0.y)**2 )**0.5
        print "________________I'm in A02"
            #if 
        self.A0 = self.A1.copy()
        self.i = self.i + self.delta_i

        return


    def find_C0(self):
     

        #self.delta_j = self.init_delta_j
        self.C1 = Node(self.CD.eval(self.j + self.delta_j),"")

        C0A0 =  ((self.C0.x - self.A0.x)**2 + (self.C0.y - self.A0.y)**2 )**0.5
        C1A0 =  ((self.C1.x - self.A0.x)**2 + (self.C1.y - self.A0.y)**2 )**0.5

        while C1A0 - C0A0 > 1e-10 :
            self.delta_j =  0.5 * self.delta_j
            self.C1 = Node(self.CD.eval(self.j + self.delta_j),"")
            print "__________________Im' here C0"
            if np.abs(self.delta_j) < self.r_t:
                #self.S_j = -1.
                self.delta_j = -self.init_delta_j + self.delta_j#-self.delta_j
            print self.delta_j
            C1A0 =  ((self.C1.x - self.A0.x)**2 + (self.C1.y - self.A0.y)**2 )**0.5
        

        self.C0 = self.C1.copy()
        self.j = self.j + self.delta_j
        return

    def get_point(self):
 
        print "__________________Im' here 3"       
        c_n = 0    # initialize the control number used to switch between different function
        while c_n < self.max_iter:
            if (c_n%2) == 0:
                self.find_A0()
            else:
                self.find_C0()
            distance = ((self.A0.x - self.C0.x)**2 + (self.A0.y - self.C0.y)**2 )**0.5
            print "distance" , distance,self.d_t

            if np.abs(distance - self.d_t) < 1.e-8:
                point = self.A0.copy()
                print "========== I'm here=======", point
                return point
            else:
                c_n = c_n + 1
            print "c_n",c_n
        if c_n == self.max_iter:
            print "Nothing happend ==============="
        return 

##
##############################################################################
##############################################################################


class offset():
    """
    class that creates a line, but offset by distance
    Inputs:
        Node0 - starting point of new line
        Node1 - end point of new line
        Path - Path object that will be stretched
    Note, currentyl only implemented in 2-D with rotation about Zaxis
    """
    def __init__(self,Node0,Node1,Path):    
        self.Node0 = Node0
        self.Node1 = Node1
        self.Path = Path
        self.os0 = ( (Node0.x-Path.eval(0.).x)**2 + (Node0.y-Path.eval(0.).y)**2 )**0.5
        self.os1 = ( (Node1.x-Path.eval(1.).x)**2 + (Node1.y-Path.eval(1.).y)**2 )**0.5 
    #
    def Path_fun(self,t):
        # get gradient
        os = (1.-t)*self.os0 + t*self.os1
        if t == 0.:
            x = self.Node0.x
            y = self.Node0.y
        elif t == 1.:
            x = self.Node1.x
            y = self.Node1.y
        else:
            # get gradient
            dt = 0.01
            if t < dt:
                t_minus = 0.    
                t_plus = 2.*t
            elif t > (1-dt):
                t_minus = 1. - 2.*t  
                t_plus = 1.
            else:
                t_minus = t - dt
                t_plus = t + dt
            plus = self.Path.eval(t_plus); minus = self.Path.eval(t_minus)
            theta = np.arctan2(plus.y-minus.y,plus.x-minus.x)
            x = self.Path.eval(t).x - os * np.sin(theta)
            y = self.Path.eval(t).y + os * np.cos(theta)
        Z = self.Node0.z # put self.Node0.z here is mainly for making the offset curve in the same plane with original curve.
        return x, y, Z
    #
    def Line(self):
        return PyFunctionPath(self.Path_fun)
##
##############################################################################
##############################################################################
class move_and_stretch():
    """
    class that creates a moved and stretched version of an existing Path element. 
    The output is a pyfunctionpath that can be used by e3prep.py
    Inputs:
        Node0 - starting point of new line
        Node1 - end point of new line
        Path - Path object that will be stretched
    Note, currentyl only implemented in 2-D with rotation about Zaxis
    """
    def __init__(self,Node0,Node1,Path):
        self.Node0 = Node0
        self.Node1 = Node1
        self.Path = Path
        self.X0 = Path.eval(0.).x
        self.Y0 = Path.eval(0.).y
        self.X1 = Path.eval(1.).x
        self.Y1 = Path.eval(1.).y 
        Length_original = ( (self.X0 - self.X1)**2 + (self.Y0 - self.Y1)**2)**0.5
        Length_new = ( (self.Node0.x - self.Node1.x)**2 + (self.Node0.y- self.Node1.y)**2)**0.5
        self.Scale = Length_new / Length_original
        self.theta_0 = np.arctan2(self.Y1-self.Y0,self.X1-self.X0)
        self.theta_new = np.arctan2(self.Node1.y-self.Node0.y,self.Node1.x-self.Node0.x)
    # 
    def Path_fun(self,t):
        x = self.Path.eval(t).x - self.X0    
        y = self.Path.eval(t).y - self.Y0
        r = x*np.cos(self.theta_0) + y*np.sin(self.theta_0)    
        s = x*np.sin(self.theta_0) - y*np.cos(self.theta_0)
        R = r*self.Scale; S = s*self.Scale
        X_new = R*np.cos(self.theta_new) + S*np.sin(self.theta_new)
        Y_new = R*np.sin(self.theta_new) - S*np.cos(self.theta_new)
        X = self.Node0.x + X_new
        Y = self.Node0.y + Y_new
        return X, Y, self.Node0.z
    #
    def Line(self):
        return PyFunctionPath(self.Path_fun)
##
##############################################################################
##############################################################################
def cart2polar(x,y,z):
    """ 
    convertes cartesian cooridinates to polar coordinates
    """
    r = (x**2+y**2)**0.5    
    theta = np.arctan2(y,x)
    return r,theta,z
##
##############################################################################
##############################################################################

def polar2cart(r,theta,z):
    """ 
    convertes polar cooridinates to cartesian coordinates
    """
    x = r * np.cos(theta)
    y = r * np.sin(theta) 
    return x,y,z
##
##############################################################################
##############################################################################
# Create outlet blocks
class NGV_outlet_chamfer():
    def __init__(self,Points,Lines0,Lines1,R_out,R_exit,bl_ratio,N,Vertical_Line="None",alpha_exit=0.,arc = 0.):
        if Vertical_Line == "None":
            self.Vertical_Line = Line(Node(0.,0.,0.), Node(0.,0.,1e-3))
        else:
            self.Vertical_Line = Vertical_Line
        self.U3 = Points[0] 
        self.B1 = Points[1]
        self.T1 = Points[2]
        self.T2 = Points[3]
        self.B6 = Points[4]
        self.D3 = Points[5]
        self.C0 = Points[6]  # YOr0
        self.C1 = Points[7]  # YOr1
        self.C2 = Points[8]  # O_chamfer
        self.C3 = Points[9]  # YOr0_BL
        self.C4 = Points[10] # YOr1_BL
        self.B9 = Points[11] # XO1_BL_2
        self.Xm = Points[12]
        self.Ym = Points[13]
        self.B7 = Points[14] # XO4_BL
        self.T3 = Points[15] # XO4

        self.Lines0 = Lines0
        self.U3B1 = Lines1[0]
        self.B1T1 = Lines1[1]
        self.T2B6 = Lines1[2]
        #self.C1B6 = Lines1[2]
        self.B6D3 = Lines1[3]
        self.center = Node(0.,0.,self.U3.z,"")

        self.bl_ratio = bl_ratio*2./3.
        self.dt_dr_exit = np.tan(alpha_exit)
        self.arc = arc
        self.r_exit = R_exit
        self.dr_exit = R_exit-R_out
        self.d_theta = - np.pi * 2. / N
        # populate data
        self.get_tangents()
        self.create_lines_points()
        self.create_patches()
        ##
    def get_tangents(self):
        l1 = 0.; l0 = 0.01
        ##
        Coords = self.Lines0[0].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[0].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.U3_dtdr = (t1-t0)/(r1-r0)

        Coords = self.Lines0[1].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[1].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.B1_dtdr = (t1-t0)/(r1-r0)


        Coords = self.Lines0[2].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[2].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.T1_dtdr = (t1-t0)/(r1-r0)

        Coords = self.Lines0[3].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[3].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.T2_dtdr = (t1-t0)/(r1-r0)

        Coords = self.Lines0[4].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[4].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.B6_dtdr = (t1-t0)/(r1-r0)

        Coords = self.Lines0[5].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[5].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.D3_dtdr = (t1-t0)/(r1-r0)

        self.B9_dtdr = self.B1_dtdr
        self.C0_dtdr = self.T2_dtdr # copy the tangent from T2 to C0

        r1,t1,z1 = cart2polar(self.C3.x,self.C3.y,self.C3.z)
        r0,t0,z0 = 0.,0.,self.C3.z
        self.C3_dtdr = (t1-t0)/(r1-r0) # copy the tangent from T2 to C1
        self.C6_dtdr = self.T2_dtdr

        self.Xm_dtdr = self.U3_dtdr
        self.Ym_dtdr = self.D3_dtdr


        ##
    def curvess_old(self,start,dt_dr,dt_dr_fraction,R,L):
        """
        function to create a curved line, which blends towards radial
        start   - starting point
        dt_dr   - gradient at start (will go towards 0)
        dt_dr_fraction - gradient of dt_dr at end. If 0, then radial line
        R       - length of the line
        L       - point at which line is being evaluated
        """
        L = 1.-L
        r,t,z = cart2polar(start.x,start.y,start.z)
        r_out = r + L * R
        z_out = z
        t_out = t + dt_dr * (L - 0.5* L**2 * (1.-dt_dr_fraction)) * (R)
        x,y,z = polar2cart(r_out,t_out,z_out)
        return x, y, z  
    def curvess(self,start,dt_dr,L):
        """
        function to create a curved line, which blends towards radial
        start   - starting point
        dt_dr   - gradient at start (will go towards 0)
        L       - point at which line is being evaluated
        """

       # print "L ===" , L
        #Use function F(x) = a * x**2 + b * x + c to describe gradient 
        arc = self.arc
        #print "self.arc, self.dt_dr_exit",self.arc, self.dt_dr_exit
        #c = dt_dr
        #b = -6. * arc + 2.*self.dt_dr_exit + 4. * dt_dr
        #a = self.dt_dr_exit - b - c
        #print dt_dr, self.dt_dr_exit
        #print "a,b,c", a,b,c
        r,t,z = cart2polar(start.x,start.y,start.z)
        dR = self.r_exit - r
        #dR = self.dr_exit
        c = dt_dr # the tangent
        a = - 5. * arc / dR**3 + 3. * self.dt_dr_exit / dR**2 + 2.0 * dt_dr / dR**2
        b = self.dt_dr_exit / dR - dt_dr / dR - a * dR

        L = 1.-L
        r,t,z = cart2polar(start.x,start.y,start.z)
        r_out = r + L * dR  #self.dr_exit  JQ 22/11/2016
        z_out = z
        X = L * dR  #self.dr_exit JQ 22/11/2016
        t_out = t + (1./3.*a*X**3 + 1./2.*b*X**2 + c*X) # * self.dr_exit

        x,y,z = polar2cart(r_out,t_out,z_out)

        return x, y, z  
        ##
    def create_lines_points(self):
        # only create axial lines, which are a helix
        test = lambda L: self.curvess(self.U3,self.U3_dtdr,L)
        self.U3L = PyFunctionPath(test)

        test = lambda L: self.curvess(self.B1,self.B1_dtdr,L)
        self.B1L = PyFunctionPath(test)

        test = lambda L: self.curvess(self.T1,self.T1_dtdr,L)
        self.T1L = PyFunctionPath(test)

        test = lambda L: self.curvess(self.T2,self.T2_dtdr,L)
        self.T2L = PyFunctionPath(test)

        test = lambda L: self.curvess(self.B6,self.B6_dtdr,L)
        self.B6L = PyFunctionPath(test)
 
        test = lambda L: self.curvess(self.D3,self.D3_dtdr,L)
        self.D3L = PyFunctionPath(test)#self.U3L.copy(); self.D3L.rotate_about_zaxis(self.d_theta)

        test = lambda L: self.curvess(self.C3,self.C0_dtdr,L)
        self.C3L = PyFunctionPath(test)

        self.C0C1 = Arc(self.C1,self.C0,self.C2)
        self.C5 = Node(self.C0C1.eval(0.5),label ="C5")
        self.C5C0 = Arc(self.C5,self.C0,self.C2)
        self.C5C1 = Arc(self.C1,self.C5,self.C2)
        self.C3C4 = Arc(self.C4,self.C3,self.C2)
        self.C6 = Node(self.C3C4.eval(0.5),label = "C6")
        self.C6C3 = Arc(self.C6,self.C3,self.C2)
        self.C6C4 = Arc(self.C6,self.C4,self.C2)
        self.C6C5 = Line(self.C6,self.C5)
        

        test = lambda L: self.curvess(self.C6,self.C6_dtdr,L)
        self.C6L = PyFunctionPath(test) 

        test = lambda L: self.curvess(self.C3,self.C3_dtdr,L)
        self.C3L = PyFunctionPath(test) 

        test = lambda L: self.curvess(self.B9,self.B9_dtdr,L)
        self.B9L = PyFunctionPath(test)


        test = lambda L: self.curvess(self.Xm,self.Xm_dtdr,L)
        self.XmL = PyFunctionPath(test)

        test = lambda L: self.curvess(self.Ym,self.Ym_dtdr,L)
        self.YmL = PyFunctionPath(test)        

        self.U3e = Node(self.U3L.eval(0.))
        self.B1e = Node(self.B1L.eval(0.))  
        self.T1e = Node(self.T1L.eval(0.))
        self.T2e = Node(self.T2L.eval(0.)) 
        self.B6e = Node(self.B6L.eval(0.))
        self.D3e = Node(self.D3L.eval(0.)) 
        self.C6e = Node(self.C6L.eval(0.))
        self.B9e = Node(self.B9L.eval(0.))
        self.C3e = Node(self.C3L.eval(0.))
        self.Xme = Node(self.XmL.eval(0.))
        self.Yme = Node(self.YmL.eval(0.))
       
        self.B9B1 = Line(self.B9,self.B1)#Polyline([self.B1L],"",(1-self.bl_ratio),1.)
        self.C5C1 = Arc(self.C5,self.C1,self.C2)

        self.B9U3 = Line(self.B9,self.U3)
        #self.B9L  = Polyline([self.B1L],"",0.,1.-self.bl_ratio)
        self.B9T1 = Line(self.B9,self.T1)
        self.C0T1 = Arc(self.C0,self.T1,Node(0.,0.,self.C0.z))
        self.C3C0 = Line(self.C3,self.C0)#Polyline([self.C0L],"",(1-self.bl_ratio),1.)
        self.C3B9 = Arc(self.C3,self.B9,Node(0.,0.,self.C3.z))

        #self.C3L =  Polyline([self.C0L],"",0.,1.-self.bl_ratio)
        self.D3C6 = Line(self.D3,self.C6)#,Node(0.,0.,self.D3.z))
        self.D3eC6 = Line(self.D3e,self.C6)
        self.B9U3e = Line(self.B9,self.U3e)
        self.U3eB9 = Line(self.U3e,self.B9)

        
        self.U3l = Node(self.U3L.eval(self.bl_ratio),"")
        self.D3l = Node(self.D3L.eval(self.bl_ratio),"")
        self.U3l_D3l = Arc(self.U3l,self.D3l,self.center)
        self.U3lB9 = Line(self.U3l,self.B9)
        self.B9U3l = Line(self.B9,self.U3l)


        self.Yml = Node(self.YmL.eval(self.bl_ratio),"")
        self.Xml = Node(self.XmL.eval(self.bl_ratio),"")
        self.Xml_Yml = Arc(self.Xml,self.Yml,self.center)
        self.XmlB9 = Line(self.Xml,self.B9)
        self.B9Xml = Line(self.B9,self.Xml)







         
        r,t,z =  cart2polar(self.Xml.x,self.Xml.y,self.Xml.z)
        r0,t0,z0 = cart2polar(self.C3.x,self.C3.y,self.C3.z)
        r0 = r 
        t0 = t0 - 0.1/(360)*2*np.pi
        x,y,z = polar2cart(r0,t0,z0)
        self.C3l = Node(x,y,z,"")
        r0 = self.r_exit
        t0 = t0 - 0.1/(360)*2*np.pi
        x,y,z = polar2cart(r0,t0,z0)
        self.C3e = Node(x,y,z,"")

        #E4 and E5 control point
        r,t,z = cart2polar(self.B9.x,self.B9.y,self.B9.z)
        r0,t0,z0 = cart2polar(self.Xml.x,self.Xml.y,self.Xml.z)
        r = r0
        x,y,z = polar2cart(r,t,z)
        self.B9l = Node(x,y,z,label = "B9l")
        r = self.r_exit
        x,y,z = polar2cart(r,t,z)
        self.B9e = Node(x,y,z,label = "B9e")





        self.U3eU3l = Polyline([self.U3L],"",0.,self.bl_ratio)
        self.C3eC3l = Line(self.C3e,self.C3l)#Polyline([self.C3L],"",0.,ratio)#
        self.C3lU3l = Arc(self.C3l,self.U3l,Node(0.,0.,self.C3l.z,""))
        self.D3lC3l = Arc(self.D3l,self.C3l,Node(0.,0.,self.D3l.z,""))
        self.C3lC3  = Line(self.C3l,self.C3)#Polyline([self.C3L],"",ratio,1.)#
        self.D3lC6  = Line(self.D3l,self.C6)

        self.U3lU3 = Polyline([self.U3L],"",self.bl_ratio,1.)
        self.D3lD3 = Polyline([self.D3L],"",self.bl_ratio,1.)
        self.D3eD3l = Polyline([self.D3L],"",0.,self.bl_ratio)


        self.XmeXml = Polyline([self.XmL],"",0.,self.bl_ratio)
        self.C3lXml = Arc(self.C3l,self.Xml,Node(0.,0.,self.C3l.z,""))
        self.YmlC3l = Arc(self.Yml,self.C3l,Node(0.,0.,self.C3l.z,""))
        


        # define YmlC6 to Bezier Curve
                #M1E3
        
        A = self.C5
        B = self.C6
        D = self.Yml
        k1 = (B.y - A.y)/(B.x - A.x)
        k2 = -1 / k1
        a = k2
        b = -1
        c = -a * B.x + B.y

        # http://stackoverflow.com/questions/3306838/algorithm-for-reflecting-a-point-across-a-line
        d = (A.x + (A.y - c )* a ) / (1 + a**2)
        C_x = 2*d - A.x
        C_y = 2*d*a - A.y + 2*c
        C = Node(C_x,C_y, A.z,'')
        self.YmlC6  = Bezier([D,C,B],"",0.,1.,1)#Line(self.Yml,self.C6)





        self.XmlXm = Polyline([self.XmL],"",self.bl_ratio,1.)
        self.YmlYm = Polyline([self.YmL],"",self.bl_ratio,1.)
        self.YmeYml = Polyline([self.YmL],"",0.,self.bl_ratio)
        #print "Yml = ",self.Yml, self.Yme, self.YmL.eval(0.)
        self.YmeC3e = Arc(self.Yme,self.C3e,Node(0.,0.,self.C3e.z,""))


        self.XmlB7 = Line(self.Xml,self.B7)
        self.B7Xml = Line(self.B7,self.Xml)
        self.B9B7 = Line(self.B9,self.B7)
        self.B9lB9 = Line(self.B9l,self.B9)
        self.B9lXml = Arc(self.B9l,self.Xml, Node(0.,0.,self.B9l.z,""))
        self.B9eXme = Arc(self.B9e,self.Xme,Node(0.,0.,self.B9e.z,""))
        self.B9eB9l = Line(self.B9e,self.B9l)
        self.B7T3 = Line(self.B7,self.T3)
        self.T1T3 = Line(self.T1,self.T3)
        self.C3lB9l = Arc(self.C3l,self.B9l,Node(0.,0.,self.B9l.z,""))
        self.C3eB9e = Arc(self.C3e,self.B9e,self.center)


        # point in middle
        self.M0 = Node(self.C0T1.eval(0.5),"")
        self.M1 = Node(self.C3B9.eval(0.5),"")
        self.Ml = Node(self.C3lB9l.eval(0.5),"")
        self.Me = Node(self.C3eB9e.eval(0.5),"")
        self.M0T1 = Arc(self.M0,self.T1,self.center)
        self.M1M0 = Bezier([self.M1,self.M0],"",0.,1.,1)
        self.M1B9 = Arc(self.M1,self.B9,self.center)

        self.MlM1 = Bezier([self.Ml,self.M1],"",0.,1.,1)
        self.MlB9l = Arc(self.Ml,self.B9l,self.center)
        self.MeMl = Line(self.Me,self.Ml)
        self.MeB9e = Arc(self.Me,self.B9e,self.center)

        self.C0M0 = Arc(self.C0,self.M0,self.center)
        self.C3M1 = Arc(self.C3,self.M1,self.center)
        self.C3lMl = Arc(self.C3l,self.Ml,self.center)
        self.C3eMe = Arc(self.C3e,self.Me,self.center)



        

    def create_patches(self):

        self.E1 = make_patch( self.B9lB9, self.M1B9, self.MlM1, self.MlB9l )
        self.E1_BL = make_patch( self.B9T1, self.M0T1, self.M1M0, self.M1B9 )
        self.E1_EX = make_patch( self.B9eB9l, self.MlB9l, self.MeMl,self.MeB9e)

        self.E2 = make_patch( self.MlM1, self.C3M1, self.C3lC3,  self.C3lMl)
        self.E2_BL = make_patch(self.M1M0, self.C0M0,self.C3C0,self.C3M1)
        self.E2_EX = make_patch(self.MeMl,self.C3lMl,self.C3eC3l,self.C3eMe)

        self.E3 = make_patch( self.C3lC3, self.C6C3, self.YmlC6, self.YmlC3l )
        self.E3_BL = make_patch( self.C3C0, self.C5C0, self.C6C5, self.C6C3 )       
        self.E3_EX = make_patch( self.C3eC3l, self.YmlC3l, self.YmeYml, self.YmeC3e )

        self.E0 = make_patch( self.XmlB7, self.B9B7, self.B9lB9, self.B9lXml )
        self.E0_BL = make_patch( self.B7T3, self.T1T3, self.B9T1, self.B9B7 )
        self.E0_EX = make_patch( self.XmeXml, self.B9lXml, self.B9eB9l, self.B9eXme )
        ##self.E2 = make_patch( self.T1L, Arc(self.C0,self.T1,Node(0.,0.,self.T1.z)), self.C0L, Arc(self.C0e,self.T1e,Node(0.,0.,self.T1e.z)) )

        #print "C0, C1, C0e, C1e", self.C0, self.C1, self.C0e, self.C1e
        #print "C1L0 and C1L1", self.C1L.eval(0), self.C1L.eval(1)
        #self.E3 = make_patch( self.C6L, self.D3C6, self.D3L, Arc(self.D3e,self.C6e,Node(0.,0.,self.C6e.z)) )
        #self.E3 = make_patch( self.T2L, self.T2B6, self.B6L, Arc(self.B6e,self.T2e,Node(0.,0.,self.T2e.z)) )   
        ##self.E4 = make_patch( self.C1L, self.C1B6, self.B6L, Arc(self.B6e,self.C1e,Node(0.,0.,self.B6e.z)) )
        #self.E4 = make_patch( self.B6L, self.B6D3, self.D3L, Arc(self.D3e,self.B6e,Node(0.,0.,self.B6e.z)) ) 
        ##self.E5 = make_patch( self.B6L, self.B6D3, self.D3L, Arc(self.D3e,self.B6e,Node(0.,0.,self.B6e.z)) )

    ##
    def Vol_E0(self,r,s,t):
        data0 =  self.E0.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E0_3D(self):
        return PyFunctionVolume(self.Vol_E0) 

    ##
    def Vol_E0_BL(self,r,s,t):
        data0 =  self.E0_BL.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E0_BL_3D(self):
        return PyFunctionVolume(self.Vol_E0_BL) 

    ##
    def Vol_E0_EX(self,r,s,t):
        data0 =  self.E0_EX.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E0_EX_3D(self):
        return PyFunctionVolume(self.Vol_E0_EX) 

    ##
    def Vol_E1_BL(self,r,s,t):
        data0 =  self.E1_BL.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E1_BL_3D(self):
        return PyFunctionVolume(self.Vol_E1_BL) 
 
    ##      
    def Vol_E1(self,r,s,t):
        data0 =  self.E1.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E1_3D(self):
        return PyFunctionVolume(self.Vol_E1)
    
    ##      
    def Vol_E1_EX(self,r,s,t):
        data0 =  self.E1_EX.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E1_EX_3D(self):
        return PyFunctionVolume(self.Vol_E1_EX)

    ##   
    def Vol_E2_BL(self,r,s,t):
        data0 =  self.E2_BL.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E2_BL_3D(self):
        return PyFunctionVolume(self.Vol_E2_BL)  
    ##   
    def Vol_E2(self,r,s,t):
        data0 =  self.E2.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E2_3D(self):
        return PyFunctionVolume(self.Vol_E2)  

    ##   
    def Vol_E2_EX(self,r,s,t):
        data0 =  self.E2_EX.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E2_EX_3D(self):
        return PyFunctionVolume(self.Vol_E2_EX)  

    ##
    def Vol_E3(self,r,s,t):
        data0 =  self.E3.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E3_3D(self):
        return PyFunctionVolume(self.Vol_E3)
    ##
    def Vol_E3_BL(self,r,s,t):
        data0 =  self.E3_BL.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E3_BL_3D(self):
        return PyFunctionVolume(self.Vol_E3_BL)
    ##
    def Vol_E3_EX(self,r,s,t):
        data0 =  self.E3_EX.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E3_EX_3D(self):
        return PyFunctionVolume(self.Vol_E3_EX)
    ##
    #def Vol_E5(self,r,s,t):
    #    data0 =  self.E5.eval(r,s)
    #    data1 = self.Vertical_Line.eval(t)
    #    return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    #def eval_E5_3D(self):
    #    return PyFunctionVolume(self.Vol_E5)    


# Create outlet blocks
class NGV_outlet():
    def __init__(self,Points,Lines0,Lines1,R_out,R_exit,Vertical_Line="None",alpha_exit=0.,arc = 0.):
        if Vertical_Line == "None":
            self.Vertical_Line = Line(Node(0.,0.,0.), Node(0.,0.,1e-3))
        else:
            self.Vertical_Line = Vertical_Line
        self.U3 = Points[0] 
        self.B1 = Points[1]
        self.T1 = Points[2]
        self.T2 = Points[3]
        self.B6 = Points[4]
        self.D3 = Points[5]
        self.Lines0 = Lines0
        self.U3B1 = Lines1[0]
        self.B1T1 = Lines1[1]
        self.T2B6 = Lines1[2]
        self.B6D3 = Lines1[3]
        self.dt_dr_exit = np.tan(alpha_exit)
        self.arc = arc
        self.dr_exit = R_exit-R_out
        # populate data
        self.get_tangents()
        self.create_lines()
        self.create_points()
        self.create_patches()
        ##
    def get_tangents(self):
        l1 = 0.; l0 = 0.02
        ##
        Coords = self.Lines0[0].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[0].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.U3_dtdr = (t1-t0)/(r1-r0)
        Coords = self.Lines0[1].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[1].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.B1_dtdr = (t1-t0)/(r1-r0)
        Coords = self.Lines0[2].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[2].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.T1_dtdr = (t1-t0)/(r1-r0)
        Coords = self.Lines0[3].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[3].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.T2_dtdr = (t1-t0)/(r1-r0)
        Coords = self.Lines0[4].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[4].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.B6_dtdr = (t1-t0)/(r1-r0)
        Coords = self.Lines0[5].eval(l0)
        r0,t0,z0 = cart2polar(Coords.x,Coords.y,Coords.z)
        Coords = self.Lines0[5].eval(l1)
        r1,t1,z1 = cart2polar(Coords.x,Coords.y,Coords.z)
        self.D3_dtdr = (t1-t0)/(r1-r0)
        ##
    def curvess_old(self,start,dt_dr,dt_dr_fraction,R,L):
        """
        function to create a curved line, which blends towards radial
        start   - starting point
        dt_dr   - gradient at start (will go towards 0)
        dt_dr_fraction - gradient of dt_dr at end. If 0, then radial line
        R       - length of the line
        L       - point at which line is being evaluated
        """
        L = 1.-L
        r,t,z = cart2polar(start.x,start.y,start.z)
        r_out = r + L * R
        z_out = z
        t_out = t + dt_dr * (L - 0.5* L**2 * (1.-dt_dr_fraction)) * (R)
        x,y,z = polar2cart(r_out,t_out,z_out)
        return x, y, z  
    def curvess(self,start,dt_dr,L):
        """
        function to create a curved line, which blends towards radial
        start   - starting point
        dt_dr   - gradient at start (will go towards 0)
        L       - point at which line is being evaluated
        """
        #Use function F(x) = a * x**2 + b * x + c to describe gradient 
        arc = self.arc
        #print self.arc, self.dt_dr_exit
        #c = dt_dr
        #b = -6. * arc + 2.*self.dt_dr_exit + 4. * dt_dr
        #a = self.dt_dr_exit - b - c
        #print dt_dr, self.dt_dr_exit
        #print "a,b,c", a,b,c
        dR = self.dr_exit
        c = dt_dr
        a = - 6. * arc / dR**3 + 3. * self.dt_dr_exit / dR**2 + 3. * dt_dr / dR**2
        b = self.dt_dr_exit / dR - dt_dr / dR - a * dR

        L = 1.-L
        r,t,z = cart2polar(start.x,start.y,start.z)
        r_out = r + L * self.dr_exit
        z_out = z
        X = L * self.dr_exit
        t_out = t + (1./3.*a*X**3 + 1./2.*b*X**2 + c*X) # * self.dr_exit

        x,y,z = polar2cart(r_out,t_out,z_out)
        return x, y, z  
        ##
    def create_lines(self):
        # only create axial lines, which are a helix
        test = lambda L: self.curvess(self.U3,self.U3_dtdr,L)
        self.U3L = PyFunctionPath(test)
        test = lambda L: self.curvess(self.B1,self.B1_dtdr,L)
        self.B1L = PyFunctionPath(test)
        test = lambda L: self.curvess(self.T1,self.T1_dtdr,L)
        self.T1L = PyFunctionPath(test)
        test = lambda L: self.curvess(self.T2,self.T2_dtdr,L)
        self.T2L = PyFunctionPath(test)
        test = lambda L: self.curvess(self.B6,self.B6_dtdr,L)
        self.B6L = PyFunctionPath(test)
        test = lambda L: self.curvess(self.D3,self.D3_dtdr,L)
        self.D3L = PyFunctionPath(test)
    ##
    def create_points(self):
        self.U3e = Node(self.U3L.eval(0.))
        self.B1e = Node(self.B1L.eval(0.))  
        self.T1e = Node(self.T1L.eval(0.))
        self.T2e = Node(self.T2L.eval(0.)) 
        self.B6e = Node(self.B6L.eval(0.))
        self.D3e = Node(self.D3L.eval(0.)) 
    def create_patches(self):
        self.E0 = make_patch( self.U3L, self.U3B1, self.B1L, Arc(self.B1e,self.U3e,Node(0.,0.,self.U3e.z)) )
        self.E1 = make_patch( self.B1L, self.B1T1, self.T1L, Arc(self.T1e,self.B1e,Node(0.,0.,self.B1e.z)) )
        self.E2 = make_patch( self.T1L, Arc(self.T2,self.T1,Node(0.,0.,self.T1.z)), self.T2L, Arc(self.T2e,self.T1e,Node(0.,0.,self.T1e.z)) )
        self.E3 = make_patch( self.T2L, self.T2B6, self.B6L, Arc(self.B6e,self.T2e,Node(0.,0.,self.T2e.z)) )     
        self.E4 = make_patch( self.B6L, self.B6D3, self.D3L, Arc(self.D3e,self.B6e,Node(0.,0.,self.B6e.z)) ) 
    ##
    def Vol_E0(self,r,s,t):
        data0 =  self.E0.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E0_3D(self):
        return PyFunctionVolume(self.Vol_E0)  
    ##      
    def Vol_E1(self,r,s,t):
        data0 =  self.E1.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E1_3D(self):
        return PyFunctionVolume(self.Vol_E1)
    ##   
    def Vol_E2(self,r,s,t):
        data0 =  self.E2.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E2_3D(self):
        return PyFunctionVolume(self.Vol_E2)  
    ##
    def Vol_E3(self,r,s,t):
        data0 =  self.E3.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E3_3D(self):
        return PyFunctionVolume(self.Vol_E3)
    ##  
    def Vol_E4(self,r,s,t):
        data0 =  self.E4.eval(r,s)
        data1 = self.Vertical_Line.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def eval_E4_3D(self):
        return PyFunctionVolume(self.Vol_E4)  



##
##############################################################################
##############################################################################
class intersection_point():
    """
    calss that return the intersection point of arbitray curves
    Jianhui Qi 29/11/2016

    """


    def __init__(self, master_Path,slave_Path):

        self.AB       = master_Path
        self.CD       = slave_Path
        self.tolerance = 1.e-9  
        self.delta_i  = 0.00001
        #self.m_r      = 0.
        #self.s_r      = 0.


    def get_tangent(self,path,ratio):
        # get the local tangent of the scanning point 

        p0 = Node(path.eval(ratio),"")
        if ratio != 1.:
            p1 =  Node(path.eval(ratio + self.delta_i),"")
        else:
            p1 =  Node(path.eval(ratio - self.delta_i),"")

        if p1.x == p0.x :
            k = "ifinite"
        else:
            k = (p1.y - p0.y) / (p1.x - p0.x)

        return k

    def get_line(self,path,ratio):

        p = Node(path.eval(ratio),"")
        k = self.get_tangent(path,ratio)
        print p,k

        if k == 0.:
            k_v = "infinite"
            A = 1.
            B = 0.
            C = -p.x

        elif k == "infinite":
            k_v = 0.
            A = 0.
            B = -1.
            C = p.y

        else:
            k_v = -1/k


            # if a point and the line tangent are given, then the line passing the 
            # point is defined as:
            # y-y0 = k(x-x0), i.e. kx - y + (-kx0 + y0) = 0. then:
            A = k_v
            B = -1
            C = -k_v*p.x + p.y
        

        line_list = [A,B,C]

        return line_list


    def p_l_distance(self,point,line_list):

        A = line_list[0]
        B = line_list[1]
        C = line_list[2]

        point_d = (A*point.x + B*point.y + C) / ((A**2 + B**2)**0.5)
                   # no abs value, that is convinent for point position determination
        return point_d

    def p_project_line(self,point,path):

        low = 0.
        high = 1.
        delta = 0.05

        p  = point
        p0 = Node(path.eval(low),"")
        p1 = Node(path.eval(low+delta),"")

        #p0_k = self.get_tangent(path,low)
        #p1_k = self.get_tangent(path,(low + delta))
        for i in range(1000):
            p0_vl = self.get_line(path,low)
            p1_vl = self.get_line(path,(low+delta))


            p0_d = self.p_l_distance(p,p0_vl)
            p1_d = self.p_l_distance(p,p1_vl)

            if low > 1. :
                ratio = 0.5
                break

            if p0_d * p1_d > (self.tolerance)**2:
                low = low + delta
            elif p0_d * p1_d < -(self.tolerance)**2:
                high = low + delta
                p0_vl = self.get_line(path,low)
                p1_vl = self.get_line(path,high)
                p0_d = self.p_l_distance(p,p0_vl)
                p1_d = self.p_l_distance(p,p1_vl)
                
                for j in range(100):
                    mid = 0.5*(low + high)
                    p2_vl = self.get_line(path,mid)
                    p2_d = self.p_l_distance(p,p2_vl)

                    if (p0_d * p2_d < -(self.tolerance)**2) and (p1_d * p2_d > (self.tolerance)**2):
                        high = mid
                        p1_vl = self.get_line(path,high)
                        p1_d = self.p_l_distance(p,p1_vl)
                        print"I'm here!!!!!!!!!+++++"

                    elif (p0_d * p2_d > (self.tolerance)**2) and (p1_d * p2_d < -(self.tolerance)**2):
                        low = mid
                        p0_vl = self.get_line(path,low)
                        p0_d = self.p_l_distance(p,p0_vl)
                        print"I'm here!!!!!!!!!-----"
                    else:
                        ratio = mid
                        break
                    print"++++++++++++++j = ",j


                print (high-low)*1.e13
                p0_vl = self.get_line(path,low)
                p1_vl = self.get_line(path,high)
                p0_d = self.p_l_distance(p,p0_vl)
                p1_d = self.p_l_distance(p,p1_vl)
                p2_vl = self.get_line(path,mid)
                p2_d = self.p_l_distance(p,p2_vl)
                print "ratio??????",high,low,mid,p1_d,p1_d,p2_d

                break
            else: 
                ratio = low + 0.5*delta
                print "ratio ==========",ratio
                break
        return ratio



    def find_point(self):

        pA = Node(self.AB.eval(0.),"")
        pC = Node(self.CD.eval(0.),"")
        pt  = Node(0.5*(pA.x+pC.x), 0.5*(pA.y+pC.y),0.5*(pA.z+pC.z),"")

        for i in range(1000):

            print "step 1"

            r_A = self.p_project_line(pt,self.AB)


            r_C = self.p_project_line(pt,self.CD)
            print "step 2",r_A,r_C
            pAt = Node(self.AB.eval(r_A),"")
            pCt = Node(self.CD.eval(r_C),"")
            pt = Node(0.5*(pAt.x+pCt.x), 0.5*(pAt.y+pCt.y),0.5*(pAt.z+pCt.z),"")

            d = ((pAt.x-pCt.x)**2 +(pAt.y-pCt.y)**2)**0.5

            if d < self.tolerance: 
                p = pt
                "it is here++++++++++++++++"
                break

            if i == 1000:
                print "!!!!!!!!!!!!"
            print "OOOOOOOOOOOOOOOOO= ",i
        print "SUCEESS!   p is", p
        return p 
            
##
##############################################################################
##############################################################################
def Arc_Curve_intersection(Arc,Curve):
    tolerance = 1.e-9
    low = 0.
    high = 1.
    Arc_start = Node(Arc.eval(0.),"")
    r,t,z = cart2polar(Arc_start.x,Arc_start.y,Arc_start.z)
    i = 0
    while i < 1000:
        mid = 0.5*(low + high)
        A = Node(Curve.eval(low),"")
        B = Node(Curve.eval(high),"")
        C = Node(Curve.eval(mid),"")
        r0,t0,z0 = cart2polar(A.x,A.y,A.z)
        r1,t1,z1 = cart2polar(B.x,B.y,B.z)
        r2,t2,z2 = cart2polar(C.x,C.y,C.z)

        if r2-r > tolerance:
            high = mid

        elif r2-r < -tolerance:
            low = mid
        else:
            curve_ratio = mid
            break
        i = i +1
    #p = Node(Curve.eval(curve_ratio),"")

    return curve_ratio

















##
##############################################################################
##############################################################################
class distance_point():
    """
    class return a point from master path which has a given distance to slave path 

    Jianhui Qi
    """

    def __init__(self,master_Path,slave_Path,distance):

        self.AB       = master_Path
        self.CD       = slave_Path
        self.d        = distance
        self.tolerace = 1.e-10

    # get the tangent of the vertical line based on p0
    # the original tangent is defined on p0 and p1
    #      |line vertical to the line p0_p1
    #      |
    #  ----o----o-----
    #      |p0  P1
    #      | 
    def get_v_tangent(self,p0,p1):
     
        k = (p1.y - p0.y) / (p1.x - p0.x)
        if k == 0:
            print "FATAL ERROR: The tangent is ", k , ", which is equal to zero", 1/0
        else:
            k_v = -1/k

        return k_v

    # get the distance from point to line
    def p_l_distance(self,point,line_variables_list):

        A = line_variables_list[0]
        B = line_variables_list[1]
        C = line_variables_list[2]

        point_d = (A*point.x + B*point.y + C) / ((A**2 + B**2)**0.5)
                   # no abs value, that is convinent for point position determination
        return point_d
    
    # return the line parameters, A B C, which are used to define a line with:
    # Ax + By + C = 0
    def p_k_line(self,line_point,line_tangent):
        
        p = line_point
        k = line_tangent
        # if a point and the line tangent are given, then the line passing the 
        # point is defined as:
        # y-y0 = k(x-x0), i.e. kx - y + (-kx0 + y0) = 0. then:
        A = k
        B = -1
        C = -k*p.x + p.y

        p_k_line_list = [A,B,C]

        return p_k_line_list


    def lines_intersection(self,line_variables_list):
        """
        return a crosecction point between given line and master path.
        if there is a crosecction point, return the evaluation ratio of the master path,
        else return -1.
        """
        
        A = line_variables_list[0]
        B = line_variables_list[1]
        C = line_variables_list[2]
        l_list = line_variables_list

        low  = 0.
        high = 1.
 
        for i in range(1000):

            mid    = 0.5*(low+high)
            p_low  = Node(self.AB.eval(low),"")
            p_high = Node(self.AB.eval(high),"")
            p_mid  = Node(self.AB.eval(mid),"")

            d_low  = self.p_l_distance(p_low,l_list)
            d_high = self.p_l_distance(p_high,l_list)
            d_mid  = self.p_l_distance(p_mid,l_list)

            if i == 0:
                if d_low * d_high > 0:
                    return -1 
                

            if np.abs(d_mid) <= self.tolerace:
                p_ratio = mid
                break

            elif d_low*d_mid > 0:
                low = mid

            elif d_low*d_mid < 0:
                high = mid

            else:
                p_ratio = mid
                break
            i += 1

        return p_ratio

    def two_p_distance(self,point0,point1):

        d = ( (point0.x - point1.x)**2 + (point0.y - point1.y)**2 )**0.5
        return d

    def find_point(self):
        start = 0.0
        delta_i = 0.001
        tangent_i = 1e-6 # use a small increase ratio to get the second point to determin the tangent
        end = start + delta_i

        for num in range(1000):

            mid = 0.5*(start + end )

            C0  = Node(self.CD.eval(start),"")
            C0_temp = Node(self.CD.eval(start + tangent_i),"")
            k_C0_v   = self.get_v_tangent(C0,C0_temp)
            C0_master_r = self.lines_intersection(self.p_k_line(C0,k_C0_v))

            C1  = Node(self.CD.eval(end), "")
            C1_temp = Node(self.CD.eval(end + tangent_i),"")
            k_C1_v   = self.get_v_tangent(C1,C1_temp)           
            C1_master_r = self.lines_intersection(self.p_k_line(C1,k_C1_v))


            Cm  = Node(self.CD.eval(mid), "")
            Cm_temp = Node(self.CD.eval(mid + tangent_i),"")
            k_Cm_v  = self.get_v_tangent(Cm,Cm_temp)
            Cm_master_r = self.lines_intersection(self.p_k_line(Cm,k_Cm_v))


            if C0_master_r + C1_master_r == -2.:
                # That means both 1st and 2nd point on master path don't pass the vertical line on slave path
                # from point 1 and point 2
                #print "comes to here1"
                start = end
                end = end + delta_i


            elif (C0_master_r == -1) and (C1_master_r != -1):
                # There is no point on master path passing the 1 vertical line, but there is some point passing 
                # the second vertical line.
                #  so this loop is to find the vertical line who can pass the most begining point on master path
                #print "comes to here2"
                start_temp = start
                end_temp = end
                
                for j in range(500):


                    mid_temp = 0.5 *(start_temp + end_temp)

                    Ct = Node(self.CD.eval(mid_temp), "")
                    Ct_temp = Node(self.CD.eval(mid_temp + tangent_i),"")
                    k_Ct_v = self.get_v_tangent(Ct,Ct_temp)
                    Ct_master_r = self.lines_intersection(self.p_k_line(Ct,k_Ct_v))

                    if Ct_master_r == -1:
                        start_temp = mid_temp


                    elif Ct_master_r > self.tolerace:
                        end_temp = mid_temp
                        print "start_temp2,end_temp,Ct_master_r",start_temp,end_temp,Ct_master_r
                    else:
                        start_temp = mid_temp
                        break
                    j = j+1
                    print "j ====",j
                    if j == 500:
                        print "The loop in find_point does not reach to converge",1/0

                start = start_temp 
                # start then become the ratio with whom the vertical line could pass the most begining point
                # on master path
                end = end + delta_i

            else:

                M0 = Node(self.AB.eval(C0_master_r),"")
                #k_C0 = -1/k_C0_v
                #C0_M0=self.p_l_distance(M0,self.p_k_line(C0,k_C0))
                C0_M0 = self.two_p_distance(M0,C0)

                Mm = Node(self.AB.eval(Cm_master_r),"")
                #k_Cm = -1/k_Cm_v
                #Cm_Mm = self.p_l_distance(Mm,self.p_k_line(Cm,k_Cm))
                Cm_Mm = self.two_p_distance(Mm,Cm)
 
                M1 = Node(self.AB.eval(C1_master_r),"")
                #k_C1 = -1/k_C1_v
                #C1_M1 = self.p_l_distance(M1,self.p_k_line(C1,k_C1))
                C1_M1 = self.two_p_distance(M1,C1)

                #print "C0_M0,Cm_Mm,C1_M1",C0_M0,Cm_Mm,C1_M1

                if (C0_M0 - self.d) * (C1_M1-self.d)  >  0:
                    start = end
                    end = end + delta_i
                    #print "I'm in 1st loop", start, end

                elif (C0_M0 - self.d) * (C1_M1-self.d)  <  0:
                    start_temp = start
                    end_temp = end
                    #print "I'm in 2nd loop"

                    #print "start_temp,end_temp",start_temp,end_temp

                    for k in range(1000):
                        mid_temp = 0.5 *(start_temp + end_temp)

                        Ct = Node(self.CD.eval(mid_temp), "")
                        Ct_temp = Node(self.CD.eval(mid_temp + tangent_i),"")
                        k_Ct_v = self.get_v_tangent(Ct,Ct_temp)
                        Ct_master_r = self.lines_intersection(self.p_k_line(Ct,k_Ct_v))
                        Mt = Node(self.AB.eval(Ct_master_r),"")
                        Ct_Mt = self.two_p_distance(Mt,Ct)
                        #print "Ct_Mt,C0_M0", Ct_Mt,C0_M0
                        #print "start_temp,end_temp",start_temp,end_temp

                        if np.abs(Ct_Mt - self.d) < 10*self.tolerace:
                            slave_ratio = mid_temp
                            #print "slave_ratio",slave_ratio
                            return Ct_master_r,slave_ratio 

                        elif (Ct_Mt - self.d) * (C0_M0-self.d) > 0:
                            start_temp = mid_temp
                        elif (Ct_Mt - self.d) * (C0_M0-self.d) < 0:
                            end_temp = mid_temp
                        else:
                            slave_ratio = mid_temp
                            #print "slave_ratio",slave_ratio
                            return Ct_master_r,slave_ratio
                        #print "start_temp,end_temp22222",start_temp,end_temp
                        k = k+1
                        if k == 1000:
                            print "loops in 3rd part of find point failed", 1/0
                else:
                    if C0_M0 == self.d:
                        master_r = C0_master_r
                        slave_ratio = start
                    else:
                        master_r = C1_master_r
                        slave_ratio = end
                    return master_r,slave_ratio


            num = num +1  
            if num == 1000:
                print "FATAL ERROR: this loops without converge",1/0

        return                  

##
##############################################################################
##############################################################################

class extruded_patch():
    def __init__(self,Patch,Path):
        self.Patch = Patch
        self.Path = Path
    def temp(self,r,s,t):
        data0 = self.Patch.eval(r,s)
        data1 = self.Path.eval(t)
        return data0.x+data1.x, data0.y+data1.y, data0.z+data1.z
    def Vol(self):
        return PyFunctionVolume(self.temp)  

