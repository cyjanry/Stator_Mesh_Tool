#!/usr/bin/env python
# This Python file uses the following encoding: utf-8
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
    def __init__(self,R_out,R_throat,R_in,R_exit,R1,R2,throat,alpha,N,throat_out,bl_ratio,C_ratio,corner_flag=0 ,r_chamfer=0.,throat_out_ratio=1., theta_0=0.,Z=0. ):
        print "Initialising supersonic_mesh()"
        self.R_out = R_out
        self.R_throat = R_throat
        self.R_in = R_in
        self.R_exit = R_exit
        self.R1 = R1
        self.R2 = R2
        self.throat = throat
        self.alpha = alpha  #alpha_trailing
        self.N = N
        self.d_theta = - np.pi * 2. / N
        self.throat_out = throat_out
        self.bl_ratio = bl_ratio       
        self.corner_flag = corner_flag
        self.r_chamfer = r_chamfer #[m] define the radius of the chamfer
        self.throat_out_ratio = throat_out_ratio
        self.theta_0= theta_0
        self.Z = Z
        self.C_ratio = C_ratio
        self.R_blade = 0.9*R_exit
        self.dt_dr_exit = np.tan(alpha)
        self.calc_fixed_nodes1()
        self.calc_flex_nodes1() 
        self.Spiral_Length(self.XR1,self.YR2)
    #
    def calc_fixed_nodes1(self):
        #------------------------------------------------# 
        # calculate main construction nodes              #
        # i.e. Nodes on outlet, inlet and throat center  #
        #------------------------------------------------#

        self.X = Node(self.R_out * np.cos(self.theta_0), self.R_out * np.sin(self.theta_0),self.Z,label="X")
        #这个地方用来计算Xt与O连线相对于OX的夹角
        angle_Xt = np.arcsin(self.R_out * np.sin(self.alpha) / self.R_throat) 
        theta_1 = self.alpha-angle_Xt

        #这个地方用来计算Xi与O连线相对于OX的夹角
        angle_Xi = np.arcsin(self.R_out * np.sin(self.alpha) / self.R_in)
        theta_2 = self.alpha-angle_Xi

        self.Xi = Node(self.R_in * np.cos(self.theta_0+theta_2), self.R_in * np.sin(self.theta_0+theta_2),self.Z,label="Xi")
        self.Xt = Node(self.R_throat * np.cos(self.theta_0+theta_1), self.R_throat * np.sin(self.theta_0+theta_1),self.Z, label="Xt")
        
        # create nodes rotated by 1 blade passage
        self.Y = self.X.copy(); self.Y.label = "Y"; self.Y.rotate_about_zaxis(self.d_theta)
        self.Yi = self.Xi.copy(); self.Yi.label = "Yi"; self.Yi.rotate_about_zaxis(self.d_theta)
        self.Yt = self.Xt.copy(); self.Yt.label = "Yt"; self.Yt.rotate_about_zaxis(self.d_theta)


        # create nodes at exit.
        angel_Xe = np.arcsin(self.R_out * np.sin(self.alpha) / self.R_exit)
        theta_3  = self.alpha-angel_Xe

        self.Xe  = Node(self.R_exit * np.cos(self.theta_0+theta_3), self.R_exit*np.sin(self.theta_0+theta_3),self.Z,label = "Xe")
        self.Ye  = self.Xe.copy(); self.Ye.label = "Ye"; self.Ye.rotate_about_zaxis(self.d_theta)

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
                  self.Xt.y - 0.5*self.throat *np.cos(self.theta_0+self.alpha), self.Z,label="Xt1")
        self.Wt1 = self.Xt1.copy(); self.Wt1.rotate_about_zaxis(-self.d_theta)
        self.Xt2 = Node(self.Xt.x - 0.5*self.throat *np.sin(self.theta_0+self.alpha), 
                  self.Xt.y + 0.5*self.throat *np.cos(self.theta_0+self.alpha),self.Z, label="Xt2")


        # Center of the chamfer, XC and YC
        rx,tx,zx = cart2polar(self.X.x, self.X.y, self.X.z)
        ry,ty,zy = cart2polar(self.Y.x, self.Y.y, self.Y.z)
     
        tXC  = self.theta_0 + self.C_ratio*2.*np.pi/self.N
        tYC  = self.theta_0 - 2*np.pi/self.N + self.C_ratio*2.*np.pi/self.N
        self.R_chamfer = self.R_out + self.r_chamfer

        xXC, yXC, zXC = polar2cart(self.R_chamfer,tXC,zx)
        xYC, yYC, zYC = polar2cart(self.R_chamfer,tYC,zy)
  

        self.XC = Node(xXC,yXC,zXC,label="XC")
        self.YC = Node(xYC,yYC,zYC,label="YC")       
    
       
        # create nodes rotated by 1 blade passage
        self.YR2 = self.XR2.copy(); self.YR2.label = "YR2"; self.YR2.rotate_about_zaxis(self.d_theta)
        self.Yt2 = self.Xt2.copy(); self.Yt2.label = "Yt2"; self.Yt2.rotate_about_zaxis(self.d_theta)
        #self.YO2 = self.XO2.copy(); self.YO2.label = "YO2"; self.YO2.rotate_about_zaxis(self.d_theta)

        # 找到倒角圆与YR2圆的公切线与倒解圆的下切点
        Q1x,Q1y,a1,b1=Common_tangent(self.YC.x,self.YC.y,self.r_chamfer,self.YR2.x,self.YR2.y,self.R2).outer_tangent_points_on_master()
        self.YCt1 = Node(Q1x,Q1y,self.Z,"YCt1") #沿主圆心到次圆心连线的方向右边的切点
        #self.Q2 = Node(a2,b2,self.Z,"Q2")

        a2,b2,Q2x,Q2y=Common_tangent(self.YC.x,self.YC.y,self.r_chamfer,self.XR1.x,self.XR1.y,self.R1).outer_tangent_points_on_master()
        self.YCt2 = Node(Q2x,Q2y,self.Z,"YCt2") #沿主圆心到次圆心连线的方向左边的切点
           

        self.XCt1 = self.YCt1.copy(); self.XCt1.label="XCt1"; self.XCt1.rotate_about_zaxis(-self.d_theta)
        self.XCt2 = self.YCt2.copy(); self.XCt2.label="XCt2"; self.XCt2.rotate_about_zaxis(-self.d_theta)

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
        self.Xt1_BL = Node((1.-self.bl_ratio)*self.Xt1.x + self.bl_ratio*self.Xt.x, 
                        (1.-self.bl_ratio)*self.Xt1.y + self.bl_ratio*self.Xt.y, 
                        self.Z, label="Xt1_BL")

        self.Wt1_BL = self.Xt1_BL.copy(); self.Wt1_BL.rotate_about_zaxis(-self.d_theta)

        self.Yt2_BL = Node((1.-self.bl_ratio)*self.Yt2.x + self.bl_ratio*self.Yt.x, 
                        (1.-self.bl_ratio)*self.Yt2.y + self.bl_ratio*self.Yt.y, 
                        self.Z, label="Yt2_BL")  
        self.Xt2_BL = self.Yt2_BL.copy(); self.Xt2_BL.label = "Xt2_BL"; self.Xt2_BL.rotate_about_zaxis(-self.d_theta)  
        
        #        r       nr
        #     o-----o----------o
        #    YC   YCt1        YCt1_BL
        #   
        #     A (Ax,Ay)   B(Bx,By)  C(Cx,Cy)
        #
        #   calculate the point YCt1_BL 
        n = 1.0
        cx,cy = mirror_point(self.YC,self.YCt1,n)
        self.YCt1_BL = Node(cx,cy, label = "YCt1_BL")

        n = 1.0
        cx,cy = mirror_point(self.YC,self.YCt2,n)
        self.YCt2_BL = Node(cx,cy, label = "YCt2_BL")    


        self.XCt1_BL = self.YCt1_BL.copy(); self.XCt1_BL.label = "XCt1_BL"; self.XCt1_BL.rotate_about_zaxis(-self.d_theta)
        self.XCt2_BL = self.YCt2_BL.copy(); self.XCt2_BL.label = "XCt2_BL"; self.XCt2_BL.rotate_about_zaxis(-self.d_theta)

        return 
    
    #接下来的部分开始定义各种各样的线条

    def Line_YCt2_Xt1(self,theta_XR1,theta_YCt2,L1=0.05,L2=0.05,L3=0.5,Node_D="None",Node_E="None"):
        t1 = theta_XR1
        t2 = theta_YCt2 #如果t2是正值，那么F点往下走，反之往上走

        A = self.Xt1.copy() ; A.label= 'A_1'
        G = self.YCt2.copy(); G.label= 'G_1'
        Length = ( (A.x - G.x)**2 + (A.y-G.y)**2 )**0.5 


        B = Node(self.XR1.x - self.R1 * np.sin(self.theta_0+self.alpha+ t1), 
                  self.XR1.y + self.R1 *np.cos(self.theta_0+self.alpha+ t1), self.Z,label="B_1")
        
        # 计算YCt2与YC连线与水平的夹角
        d_y = self.YCt2.y - self.YC.y
        alpha = np.arccos(d_y/self.r_chamfer)

        F = Node(self.YC.x - self.r_chamfer * np.sin(alpha+ t2), 
                  self.YC.y + self.r_chamfer *np.cos(alpha+ t2), self.Z,label="F_1")
        # 这几个点需要重点讨论一下，和狗哥讨论讨论
        C = Node(B.x - L1*Length * np.cos(self.theta_0+self.alpha+t1),
                B.y - L1*Length * np.sin(self.theta_0+self.alpha+t1),self.Z,label="C_1")

        E = Node(F.x + L2*Length * np.cos(alpha+t2),
                F.y + L2*Length * np.sin(alpha+t2),self.Z,label="E_1")

        #print "------",(A.y - G.y)/(A.x - G.x),
        theta_AG = np.arctan2((A.y - G.y),(A.x - G.x))


        D = Node( 0.5*(A.x + G.x) - np.cos(0.5*np.pi-theta_AG)*L3*Length, 
            0.5*(A.y + G.y) - np.sin(0.5*np.pi-theta_AG)*L3*Length, self.Z, label="D_1")
        # 这个地方需要严格讨论是否要再加一个控制点，加的话，怎么加???

        if t1 != 0:
            Line2 = Arc(B,A,self.XR1)
            Line1 = Bezier([F,E,D,C,B],"",0.,1.,1)
            Line0 = Arc(G,F,self.YC)
            temp  = Polyline([Line0,Line1,Line2],"",0.,1.,1) 
        else:
            Line1 = Bezier([F,E,D,C,B],"",0.,1.,1)
            Line0 = Arc(G,F,self.YC)
            temp  = Polyline([Line0,Line1],"",0.,1.,1) 

        #Line3 = Line(self.YCt2, self.Xt1)

        self.Line_YCt2_Xt1 = temp
        self.Line_YCt2BL_Xt1BL = move_and_stretch(self.YCt2_BL,self.Xt1_BL,self.Line_YCt2_Xt1).Line()



#        A_ = A.copy(); A_.rotate_about_zaxis(self.d_theta)
#        G_ = self.XCt2.copy()
#        Length_ = ( (A_.x - G_.x)**2 + (A_.y-G_.y)**2 )**0.5 
#        B_ = B.copy(); B_.rotate_about_zaxis(self.d_theta)
#        F_ = F.copy(); F_.rotate_about_zaxis(self.d_theta)
#        C_ = C.copy(); C_.rotate_about_zaxis(self.d_theta)
#        E_ = E.copy(); E_.rotate_about_zaxis(self.d_theta)


        return 



    def Line_YCt1_Yt2(self,theta_YR2,theta_YCt1,L1=0.1,L2=0.1,Node_D="None",Node_E="None"):
        t1 = theta_YR2
        t2 = theta_YCt1

        A = self.Yt2.copy(); A.label = "A_2"
        G = self.YCt1.copy(); G.label = "G_2"

        Length = ( (A.x - G.x)**2 + (A.y-G.y)**2 )**0.5 

        

        BX = Node(self.XR2.x + self.R2 * np.sin(self.theta_0+self.alpha+t1), 
                  self.XR2.y - self.R2 *np.cos(self.theta_0+self.alpha+t1), self.Z,label="")
        B = BX.copy(); B.label="B_2"; B.rotate_about_zaxis(self.d_theta)

        # 计算YCt2与YC连线与水平的夹角
        d_y = abs(self.YCt1.y - self.YC.y)
        alpha = np.arccos(d_y/self.r_chamfer)#; print alpha,"XX"


        F = Node(self.YC.x + self.r_chamfer * np.sin(alpha- t2), 
                  self.YC.y - self.r_chamfer *np.cos(alpha- t2), self.Z,label="F_2")

        #使用XR2的位置计算CX，然后再把它转d_theta得到下面的C
        CX = Node(BX.x - L1*Length * np.cos(self.theta_0+self.alpha+t1),
                BX.y - L1*Length * np.sin(self.theta_0+self.alpha+t1),self.Z,label = "")
        C = CX.copy(); C.label = "C_2"; C.rotate_about_zaxis(self.d_theta)

        E = Node(F.x + L2*Length * np.cos(alpha+t2),
                F.y + L2*Length * np.sin(alpha+t2),self.Z,label="E_2")

        #combine lines for output

        if t1 != 0.:
                Line2 = Arc(B,A,self.YR2)
                Line1 = Bezier([F,E,C,B],"",0.,1.,1)
                Line0 = Arc(G,F,self.YC)
                temp  = Polyline([Line0,Line1,Line2],"",0.,1.,1)
        else:
                Line1 = Bezier([F,E,C,B],"",0.,1.,1)
                Line0 = Arc(G,F,self.YC)
                temp  = Polyline([Line0,Line1],"",0.,1.,1)            

        self.Line_YCt1_Yt2 = temp # TODO: 这一块也要好好考虑一下

        self.Line_YCt1BL_Yt2BL = move_and_stretch(self.YCt1_BL,self.Yt2_BL,temp).Line()

        A_ = A.copy(); A_.rotate_about_zaxis(-self.d_theta)
        G_ = self.XCt2.copy()
        Length_ = ( (A_.x - G_.x)**2 + (A_.y-G_.y)**2 )**0.5 
        B_ = B.copy(); B_.rotate_about_zaxis(-self.d_theta)
        F_ = F.copy(); F_.rotate_about_zaxis(-self.d_theta)
        C_ = C.copy(); C_.rotate_about_zaxis(-self.d_theta)
        E_ = E.copy(); E_.rotate_about_zaxis(-self.d_theta)        


        self.Line_XCt1_Xt2 = temp.copy(); self.Line_XCt1_Xt2.rotate_about_zaxis(-self.d_theta)
        self.Line_XCt1BL_Xt2BL = move_and_stretch(self.XCt1_BL, self.Xt2_BL, self.Line_XCt1_Xt2).Line()


        return 



    def calc_depend_nodes(self):
        # function that calcualte the nodes that depend on other curves.
        # 些函数是用来计算依赖于其它直线的点。
        # 这个函数计算出来的点不影响geometry的设计，但是会影响网格质量的好与坏
        # 总得来说，要想使用这个函数得到点，要首先得到依存的直线

        length_between_YCt1_Yt2 = ( (self.YCt1.x - self.Yt2.x)**2 + (self.YCt1.y - self.Yt2.y)**2 )**0.5
        length_between_YCt2_Xt1 = ( (self.YCt2.x - self.Xt1.x)**2 + (self.YCt2.y - self.Xt1.y)**2 )**0.5
        #length_between_XCt1_Xt2 = ( (self.XCt1.x - self.Xt2.x)**2 + (self.XCt1.y - self.Xt2.y)**2 )**0.5

        self.l_ratio_Xm = length_between_YCt1_Yt2/length_between_YCt2_Xt1

        # 因为是反方向的，所以要相减
        self.Xm = Node(self.Line_YCt2_Xt1.eval(1.-self.l_ratio_Xm),label="Xm")
        self.Xm_BL = Node(self.Line_YCt2BL_Xt1BL.eval(1.-self.l_ratio_Xm),label="Xm_BL")

        # 接下来找到XC与Xm连线(Beizer_Curve)与X_Xi的交点。
        self.Line_XCt1_XmBL = Line(self.XCt1_BL,self.Xm_BL)
        self.Line_X_Xi  = Line(self.X,self.Xi)

        TEST = distance_point(self.Line_X_Xi,self.Line_XCt1_XmBL,1e-5)
        line_list = TEST.line_list_two_points(self.XCt1,self.Xm_BL)
        #print line_list
        master_ratio = TEST.lines_intersection(line_list)
        #print master_ratio
        self.XB = Node(self.Line_X_Xi.eval(master_ratio),label = "XB")
        self.YB = self.XB.copy(); self.YB.label = "YB"; self.YB.rotate_about_zaxis(self.d_theta)


        Arc_YCt2_YCt1 = Arc(self.YCt2,self.YCt1,self.YC)
        self.YCt0 = Node(Arc_YCt2_YCt1.eval(0.5),label = "YCt0")

        Arc_YCt2BL_YCt1BL = Arc(self.YCt2_BL,self.YCt1_BL,self.YC)
        self.YCt0_BL = Node(Arc_YCt2BL_YCt1BL.eval(0.5),label = "YCt0_BL")

        self.XCt0 = self.YCt0.copy(); self.XCt0.label = "XCt0"; self.XCt0.rotate_about_zaxis(-self.d_theta)
        self.XCt0_BL = self.YCt0_BL.copy(); self.XCt0_BL.label = "XCt0_BL"; self.XCt0_BL.rotate_about_zaxis(-self.d_theta)




        #在这儿使用出口曲线函数计算XeL的位置
        r0,t0,z0 = cart2polar(self.X.x,self.X.y,self.X.z)
        r1,t1,z1 = cart2polar(self.Xt.x,self.Xt.y,self.Xt.z)
        X_dtdr = (t1-t0)/(r1-r0)
        arc1 = -3./180.*np.pi        
        test = lambda L: self.curvess(self.X,X_dtdr,arc1,L)
        temp = PyFunctionPath(test)
        self.Line_XeL_X = Polyline([temp],"",0,1.,1) 
        self.XeL = Node(self.Line_XeL_X.eval(0.),label = "XeL")


        r0,t0,z0 = cart2polar(self.Y.x,self.Y.y,self.Y.z)
        r1,t1,z1 = cart2polar(self.Yt.x,self.Yt.y,self.Yt.z)
        Y_dtdr = (t1-t0)/(r1-r0)
        arc2 = -3./180.*np.pi
        test = lambda L: self.curvess(self.Y,Y_dtdr,arc2,L)
        temp = PyFunctionPath(test)        
        self.Line_YeL_Y = Polyline([temp],"",0,1.,1) 
        self.YeL = self.XeL.copy(); self.YeL.label = "YeL"; self.YeL.rotate_about_zaxis(self.d_theta)



        # 使用XeL_Xt的长度来估算Xn和Xn_BL的位置
        length_between_XeL_Xt =  ( (self.XeL.x - self.Xt.x)**2 + (self.XeL.y - self.Xt.y)**2 )**0.5
        self.l_ratio_Xn = length_between_XeL_Xt/length_between_YCt2_Xt1

        self.Xn = Node(self.Line_YCt2_Xt1.eval(1.-self.l_ratio_Xn),label = "Xn")
        self.Xn_BL = Node(self.Line_YCt2BL_Xt1BL.eval(1.-self.l_ratio_Xn),label = "Xn_BL")




        #使用极坐标转换得到 YCt1_YCt1BL 出口延长线在exit圆弧上的切点
        r,t,z = cart2polar(self.YCt2_BL.x,self.YCt2_BL.y,self.Z)
        r = self.R_exit; x,y,z = polar2cart(r,t,z)
        self.YCte = Node(x,y,z,label = "YCte")
        self.XCte = self.YCte.copy(); self.XCte.rotate_about_zaxis(-self.d_theta);self.XCte.label = "XCte"
        #print "YCte's r=",r

        # control points at other sides of YCt1_BL        
        #        r       nr
        #     o-----o----------o
        #   YCt2  YCt2_BL  YCt2_ctrl
        #   
        #     A (Ax,Ay)   B(Bx,By)  C(Cx,Cy)
        #
        #   calculate the point YCt1_BL 
        n = 1.0
        cx,cy = mirror_point(self.YCt2,self.YCt2_BL,n)
        self.YCt2_ctrl = Node(cx,cy,self.Z, label = "YCt2_ctrl")

        #做出口的弯曲线

        r0,t0,z0 = cart2polar(self.YCt0.x,self.YCt0.y,self.YCt0.z)
        r1,t1,z1 = cart2polar(self.YCt0_BL.x,self.YCt0_BL.y,self.YCt0_BL.z)
        YCt0_BL_dtdr = (t1-t0)/(r1-r0)
        arc3 = -1/180.*np.pi
        test = lambda L: self.curvess(self.YCt0_BL,YCt0_BL_dtdr,arc3,L)
        temp = PyFunctionPath(test)
        self.Line_YCt0BL_YCt0e = Polyline([temp],"",1.0,0.,1) #将这条线反过来
        self.YCt0e = Node(self.Line_YCt0BL_YCt0e.eval(1.),label = "YCt0e")
        self.XCt0e = self.YCt0e.copy(); self.XCt0e.rotate_about_zaxis(-self.d_theta);self.XCt0e.label = 'XCt0e'


        n = 1.0
        cx,cy = mirror_point(self.YCt0,self.YCt0_BL,n)
        self.YCt0_ctrl = Node(cx,cy,self.Z, label = "YCt0_ctrl")

        # 使用Beizer Curve找到直线X Xi上的Xa点
        
        n = 2.0
        cx,cy = mirror_point(self.YCt0,self.YCt0_BL,n)
        self.Ya2 = Node(cx,cy,self.Z, label = "Ya2")

        self.XCt0_ctrl = self.YCt0_ctrl.copy(); self.XCt0_ctrl.label = "XCt0_ctrl"; self.XCt0_ctrl.rotate_about_zaxis(-self.d_theta)
        self.Xa2 = self.Ya2.copy(); self.Xa2.label = "Xa2"; self.Xa2.rotate_about_zaxis(-self.d_theta)




        k = (self.Xt.y - self.X.y)/ (self.Xt.x - self.X.x)
        kv = -1./k
        temp_l = Line(self.Xe,self.Xt)
        TEST = distance_point(temp_l,temp_l,1e-6)
        list_l = TEST.p_k_line(self.Xa2,kv)
        master_ratio = TEST.lines_intersection(list_l)

        self.Xa = Node(temp_l.eval(master_ratio),label="Xa")
        self.Ya = self.Xa.copy(); self.Ya.label = "Ya"; self.Ya.rotate_about_zaxis(self.d_theta)


        #然后用Xt_Xa的长度来和Xt1_Yc的长度进行对比，得到Xa1和Xa1_BL的位置
        length_between_Xa_Xt = ( (self.Xa.x - self.Xt.x)**2 + (self.Xa.y - self.Xt.y)**2 )**0.5
        self.l_ratio_Xa1          = length_between_Xa_Xt / length_between_YCt2_Xt1
        self.Xa1 = Node(self.Line_YCt2_Xt1.eval(1. -self.l_ratio_Xa1), label = "Xa1")
        self.Xa1_BL = Node(self.Line_YCt2BL_Xt1BL.eval(1.-self.l_ratio_Xa1), label = "Xa1_BL") 

        #print self.Line_XeL_X.eval(0.), self.Line_XeL_X.eval(1.),self.X


        #在这儿加上出口的一些点

        r,t,z  = cart2polar(self.XCte.x, self.XCte.y, self.Z)
        x,y,z  = polar2cart(self.R_blade,t,z)
        self.XCte_O = Node(x,y,z,label = "XCte_O")

        r,t,z  = cart2polar(self.XCt0e.x,self.XCt0e.y, self.Z)
        x,y,z  = polar2cart(self.R_blade,t,z)
        self.XCt0e_O = Node(x,y,z,"XCt0e_O")

        r,t,z  = cart2polar(self.XeL.x, self.XeL.y, self.Z)
        x,y,z  = polar2cart(self.R_blade,t,z)
        self.XeL_O  = Node(x,y,z, "XeL_O")

        self.YCte_O = self.XCte_O.copy(); self.YCte_O.rotate_about_zaxis(self.d_theta); self.YCte_O.label = "YCte_O"
        self.YCt0e_O= self.XCt0e_O.copy(); self.YCt0e_O.rotate_about_zaxis(self.d_theta); self.YCt0e_O.label = "YCt0e_O"
        self.YeL_O = self.XeL_O.copy(); self.YeL_O.rotate_about_zaxis(self.d_theta); self.YeL_O.label = "YeL_O"




        return

   #


    def Line_XR1t_YR2t(self,X1s_ratio,TOP_ratio,theta_XT1 = -1.6,theta_YR2 = 0.4):
        
        t1 = theta_XT1
        t2 = theta_YR2

        #TOP_ratio = 0.25
        #print "-------------------------------"
        #print "Initializing curve XR1t_YR2t..."
        #print "-------------------------------"

        nose_path_temp = PyFunctionPath(self.Path_Xt1_Yt2_3)



        A = self.Xt1.copy();A.label = "A_3"
        L = self.Yt2.copy();L.label = "L_3"
        #B = nose_path_temp.eval(0.1)
        B = Node(self.XR1.x - self.R1 * np.sin(self.theta_0+self.alpha + t1), 
                  self.XR1.y + self.R1 *np.cos(self.theta_0+self.alpha + t1), self.Z,label="B_3")


        temp1 = Arc(A,B,self.XR1)


        d_y = B.y - self.XR1.y
        alpha = np.arccos(d_y/self.R1)
        Length =   ( (B.x - L.x)**2 + (B.y - L.y)**2 )**0.5

        C = Node(B.x + 0.05*Length * np.cos(alpha ) , 
                  B.y-0.05*Length * np.sin(alpha) , self.Z,label="C_3")


        #C = nose_path_temp.eval(0.25);C.label='c'
        D = Node(nose_path_temp.eval(0.3));D.label = 'D_3'
        E = Node(nose_path_temp.eval(0.4));E.lable = 'E_3'
        F = Node(nose_path_temp.eval(0.45));F.label = 'F_3'
        G = Node(nose_path_temp.eval(0.5));G.label = 'G_3'
        H = Node(nose_path_temp.eval(0.55));H.label = 'H_3'

        
        KX = Node(self.XR2.x + self.R1 * np.sin(self.theta_0+self.alpha+t2), 
                  self.XR2.y - self.R1 *np.cos(self.theta_0+self.alpha+t2), self.Z,label="")
        K = KX.copy(); K.label="K_3"; K.rotate_about_zaxis(self.d_theta)
        

        theta_HK = np.arctan2((H.y - K.y),(H.x - K.x))      

        I = Node(0.5*(H.x +K.x)- np.cos(0.5*np.pi-theta_HK)*0.1*Length,
                 0.5*(H.y + K.y)- np.sin(0.5*np.pi-theta_HK)*0.1*Length, self.Z, label = 'I_3')


        #D = Node( 0.5*(A.x + G.x) - np.cos(0.5*np.pi-theta_AG)*0.1*Length, 
        #    0.5*(A.y + G.y) - np.sin(0.5*np.pi-theta_AG)*0.1*Length, self.Z, label="D_1")


        # 计算YCt2与YC连线与水平的夹角
        d_y = K.y - self.YR2.y
        alpha = np.arccos(d_y/self.R1)



        J = Node(K.x - 0.1*Length * np.cos(alpha),
                K.y + 0.1*Length * np.sin(alpha),self.Z,label="J_3") 



        temp2 = Arc(K,L,self.YR2)

        self.nose_path = Polyline([temp1, Bezier([B,C,D,E,F,G,H,I,J,K],""''"",0.,1.,1), temp2] , "",0. , 1., 1)



        self.nose_path_bl= offset(self.Xt1_BL , self.Yt2_BL , self.nose_path).Line()
        # The following code is used to calculate XR1s, XR1s_BL0 and XR1s_BL
        self.rotated_nose_path = self.nose_path.copy().rotate_about_zaxis(-self.d_theta)# PyFunctionPath(self.Path_Wt1_Xt2_3)
        self.rotated_nose_path_bl = offset(self.Wt1_BL , self.Xt2_BL, self.rotated_nose_path).Line()


        self.XR1s     =  Node(self.nose_path.eval(X1s_ratio), label= "XR1s")
        self.WR1s     =  self.XR1s.copy(); self.WR1s.rotate_about_zaxis(-self.d_theta)

        self.TOP      =  Node(self.nose_path.eval(TOP_ratio), label= "TOP")
        self.WTOP     =  Node(self.rotated_nose_path.eval(TOP_ratio), label= "WTOP")#self.TOP.copy(); self.WTOP.rotate_about_zaxis(-self.d_theta);self.WTOP.label = "WTOP"

        self.XR1s_BL0 =  Node(self.nose_path_bl.eval(X1s_ratio),label = "XR1s_BL0")
        self.TOP_BL   =  Node(self.nose_path_bl.eval(TOP_ratio),label = "TOP_BL")
        self.WTOP_BL  =  Node(self.rotated_nose_path_bl.eval(TOP_ratio),label = "WTOP_BL") 
        #self.WTOP_BL = self.TOP_BL.copy(); self.WTOP_BL.rotate_about_zaxis(-self.d_theta)


        self.L1_BL_S   = Polyline([self.nose_path],"line1",0.,X1s_ratio)
        self.TOPL_BL_W = Polyline([self.nose_path],"line2",TOP_ratio,X1s_ratio)

        #        d       nd
        #     o-----o----------o
        #     A     B          C
        #   
        #     A (Ax,Ay)   B(Bx,By)  C(Cx,Cy)
        #
        #   calculate the point XR1s_BL         
        n = 2.0   # set how far the point c should located.
        cx, cy = mirror_point(self.XR1s,self.XR1s_BL0,n)
        self.XR1s_BL = Node(cx,cy, label = "XR1s_BL")
        self.WR1s_BL = self.XR1s_BL.copy(); self.WR1s_BL.rotate_about_zaxis(-self.d_theta)

        # Then calculation line Xt1BL_XR1sBL
        self.Xt1BL_XR1sBL = move_and_stretch(self.Xt1_BL,self.XR1s_BL,self.L1_BL_S).Line() 
        self.TOPBL_XR1sBL = move_and_stretch(self.TOP_BL,self.XR1s_BL,self.TOPL_BL_W).Line()
        self.WTOPBL_WR1sBL =  move_and_stretch(self.WTOP_BL,self.WR1s_BL,self.TOPL_BL_W).Line()
        
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
    def spiral_angle(self,r_start,r_end,t_start,t_end,t,dt=0.0001):
        # function to get angle of sprial line
        t_p = t+dt
        if t_p > 1.: 
            t_p =1.
        t_m = t-dt
        if t_m < 0.:
            t_m = 0.
        #print t_m
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
        start = self.XR1.copy();start.label="" # center of circle at start
        end = self.YR2.copy();end.label="" # center of circle at end        #

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
        #print "------?----------"
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
        #print "XX",L0
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
        #print "Inseide function:",t,tt, x,y,z
        return (x,y,z)

    def Path_Wt1_Xt2_3(self,t):
        #PyFunctionPath(t) to evaluate line
        WR1 = self.XR1.copy(); WR1.rotate_about_zaxis(-self.d_theta)
        Wt1 = self.Xt1.copy(); Wt1.rotate_about_zaxis(-self.d_theta)
        start = WR1.copy();start.label="" # center of circle at start
        end = self.XR2.copy();end.label="" # center of circle at end        #

        r_start = np.sqrt(start.x**2 + start.y**2)
        t_start = np.arctan2(start.y,start.x)
        r_end = np.sqrt(end.x**2 + end.y**2)
        t_end = np.arctan2(end.y,end.x)
        #
        R_start = ( (Wt1.y - start.y)**2 + (Wt1.x - start.x)**2 )**0.5
        R_end   = ( (self.Xt2.y - end.y  )**2 + (self.Xt2.x - end.x  )**2 )**0.5
        #
        theta_start = np.arctan2(Wt1.y - start.y ,Wt1.x - start.x)
        theta_end   = np.arctan2(self.Xt2.y - end.y  ,self.Xt2.x - end.x  )

        theta_mstart = self.spiral_angle(r_start,r_end,t_start,t_end,0.)+np.pi/2.
        theta_mend = self.spiral_angle(r_start,r_end,t_start,t_end,1.)+np.pi/2.

        if theta_start < theta_mstart:
            theta_start += 2.*np.pi
        #print "Angles after:", theta_start, theta_mstart, theta_mend, theta_end
        L0 = R_start * (theta_start-theta_mstart)

        L1 = self.L1
        #print "XX",L0
        if theta_mend >= theta_end:
            L2 = R_end * (theta_mend - theta_end)
            LL = L0+L1+L2
            #print L0, L1, L2, LL
            if t < L0/LL:
                tt = t * LL/L0
                C_x = start.x 
                C_y = start.y 

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
        #print "Inseide function:",t,tt, x,y,z
        return (x,y,z)



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




    def Out_let_nodes(self):
        #这一部分用来定义出口的点
    #    self.XO = 
        return

    def curvess(self,start,dt_dr,arc,L):
        """
        function to create a curved line, which blends towards radial
        start   - starting point
        dt_dr   - gradient at start (will go towards 0)
        L       - point at which line is being evaluated
        """
        #print "L:", L
        #Use function F(x) = a * x**2 + b * x + c to describe gradient 
        #arc = -3./180.*np.pi#self.arc

        #print "self.arc, self.dt_dr_exit",self.arc, self.dt_dr_exit
        #c = dt_dr
        #b = -6. * arc + 2.*self.dt_dr_exit + 4. * dt_dr
        #a = self.dt_dr_exit - b - c
        #print dt_dr, self.dt_dr_exit
        #print "a,b,c", a,b,c
        r,t,z = cart2polar(start.x,start.y,start.z)

        dR = self.R_exit - r
        #dR = self.dr_exit

        c = dt_dr # the tangent
        a = - 6. * arc / dR**3 + 5. * self.dt_dr_exit / dR**2 + 2.0 * dt_dr / dR**2
        b = self.dt_dr_exit / dR - dt_dr / dR - a * dR

        L = 1.-L
        r,t,z = cart2polar(start.x,start.y,start.z)
        r_out = r + L * dR  #self.dr_exit  JQ 22/11/2016
        z_out = z
        X = L * dR  #self.dr_exit JQ 22/11/2016
        t_out = t + (1./3.*a*X**3 + 1./2.*b*X**2 + c*X) # * self.dr_exit

        x,y,z = polar2cart(r_out,t_out,z_out)
        #print x,y,z

        return x, y, z 





##
##############################################################################
##############################################################################
def mirror_point(point_A, point_B, n):
    # A, B, d, n is known. C should be determined.
    #
    #        d       nd
    #     o-----o----------o
    #     A     B          C
    #   
    #     A (Ax,Ay)   B(Bx,By)  C(Cx,Cy)
    #
    #   calculate the point XR1s_BL 

    A = point_A
    B = point_B
    
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

    return cx,cy



    
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
            dt = 0.001 # This is important
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
class Common_tangent():
    # This class is used to calculate the common tangents for two circles.
    # 这个类是用来计算两个给定圆的公切线，及切点

    def __init__(self,x1,y1,r1,x2,y2,r2):
        self.x1  = x1
        self.y1  = y1
        self.r1  = r1
        self.x2  = x2
        self.y2  = y2
        self.r2  = r2
        self.outer_tangent()

    def outer_tangent(self):
        #这个函数的功能是返回两个给定圆的公切线的直线参数 k,b.  y = kx+b
        # x1,y1 是主圆的圆心坐标 r1 是主圆的半径
        x1  = self.x1
        y1  = self.y1
        r1  = self.r1
        x2  = self.x2
        y2  = self.y2
        r2  = self.r2


        #http://bbs.emath.ac.cn/thread-5803-1-1.html
        kA = (-x2+x1+r1-r2) * (x2-x1+r1-r2)
        kB = 2*(y1-y2) * (x1-x2)
        kC = (y2+r1-y1-r2) * (-y2+r1+y1-r2)

        k1 = (-kB + np.sqrt(kB*kB - 4.*kA*kC )) / (2.*kA)
        k2 = (-kB - np.sqrt(kB*kB - 4.*kA*kC )) / (2.*kA)

        bA = (-x2+x1+r1-r2)*(x2-x1+r1-r2)
        bB = (2*r1*r1*y2 - 2*r1*r2*y1 - 2*r1*r2*y2 + 2*r2*r2*y1 - 2*x1*x1*y2+2*x1*x2*y1+2*x1*x2*y2 - 2*x2*x2*y1)
        bC = r1*r1*x2*x2+r1*r1*y2*y2 - 2*r1*r2*x1*x2 - 2*r1*r2*y1*y2 + r2*r2*x1*x1 + r2*r2*y1*y1 - x1*x1*y2*y2+ 2*x1*x2*y1*y2 - x2*x2*y1*y1
        
        b1 = (-bB + np.sqrt(bB*bB - 4.*bA*bC)) / (2.*bA)
        b2 = (-bB - np.sqrt(bB*bB - 4.*bA*bC)) / (2.*bA)
        self.outer_tangent_k1 = k1
        self.outer_tangent_k2 = k2
        self.outer_tangent_b1 = b1
        self.outer_tangent_b2 = b2

        # 同大和同小组合，给定两条直线的方程
        """
        if k1 > k2:
            self.outer_tangent_k1 = k1
            self.outer_tangent_k2 = k2
        else:
            self.outer_tangent_k1 = k2
            self.outer_tangent_k2 = k1

        if b1 > b2:
            self.outer_tangent_b1 = b1
            self.outer_tangent_b2 = b2
        else:
            self.outer_tangent_b1 = b2
            self.outer_tangent_b2 = b1

        """
        #print k1, b1, k2, b2
        return 



    def outer_tangent_points_on_master(self):

        #这个函数的功能是，找出两条外公切线与主圆的两个切点
        #沿着主圆心到次圆心的方向，a1b1是在右手边的切点， a2b2是在左手边的切点。所以要好好使用这个特性
        # a1,b1 

        a1 =    (self.x1 + (self.outer_tangent_b1 + self.y1)*self.outer_tangent_k1) /(1 + self.outer_tangent_k1**2)
        b1 =    self.outer_tangent_k1*a1 - self.outer_tangent_b1

        a2 =    (self.x1 + (self.outer_tangent_b2 + self.y1)*self.outer_tangent_k2) /(1 + self.outer_tangent_k2**2)
        b2 =    self.outer_tangent_k2*a2 - self.outer_tangent_b2
        return a1, b1, a2, b2

##
##############################################################################
##############################################################################




##
##############################################################################
##############################################################################
class distance_point():
    
    #class return a point from master path which has a given distance to slave path 
    #这个类是用来找到一个点。 首先给定两条线，主线和次线。 这个类在主线上找到一个点，使这个点到次线的
    #距离为给定值。
    #Jianhui Qi/ 齐建荟
    

    def __init__(self,master_Path,slave_Path,distance):

        self.AB       = master_Path  # 定义主线
        self.CD       = slave_Path   # 定义次线
        self.d        = distance     # 给定距离
        self.tolerace = 1.e-10       # 给定收敛公差

    # get the tangent of the vertical line based on p0
    # the original tangent is defined on p0 and p1
    #      |line vertical to the line p0_p1
    #      |
    #  ----o----o-----  <-- k
    #      |p0  P1
    #      | 
    # 
    #      ^
    #      :
    #     k_v



    def get_v_tangent(self,p0,p1):
        # 这个函数用来返回垂直于p0 p1连线的直线的斜率
     
        k = (p1.y - p0.y) / (p1.x - p0.x) # 这一步是找到line_p0_p1的斜率
        if k == 0:
            print "FATAL ERROR: The tangent is ", k , ", which is equal to zero", 1/0
        else:
            k_v = -1/k   # 给定line_p0_p1的垂线斜率

        return k_v

    def line_list_two_points(self,p0,p1):
        k = (p1.y - p0.y) / (p1.x - p0.x) # 这一步是找到line_p0_p1的斜率
        if k == 0:
            print "FATAL ERROR: The tangent is ", k , ", which is equal to zero", 1/0
        else:
            k_v = -1/k   # 给定line_p0_p1的垂线斜率 

        b = p0.y - k*p0.x
        return [k, -1, b]

    def p_l_distance(self,p,line_variables_list):
        # get the distance from point to line
        # 此函数用来得到给定点到给定直线的距离
        # 其中直接议程为  Ax + By + C = 0, 所以直接给定 [A, B, C]这样一个直线参数列表
        # 给出来的距离并不是距离的绝对值，这样有助于决定点和直线的相对位置。
        A = line_variables_list[0]
        B = line_variables_list[1]
        C = line_variables_list[2]

        point_d = (A*p.x + B*p.y + C) / ((A**2 + B**2)**0.5)
                   # no abs value, that is convinent for point position determination
        return point_d
    

    def p_k_line(self,line_point,line_tangent):
        # 此函数用来得到一条直线方程的参数。
        # 这条直线是给定点和斜率之后，转换为 Ax + By + C = 0的形式
        # return the line parameters, A B C, which are used to define a line with:
        # Ax + By + C = 0

        
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
        #此函数用来找到一条直线与主线（可以是直线，弧线，贝赛尔曲线）的交点在主线上的位置比例
        #首先给定次线的直线方程，是列表的形式。
        #主线已经在类的开使定义了，为self.AB， 然后通过二分法求解交点
        #return a crosecction point between given line and master path.
        #if there is a crosecction point, return the evaluation ratio of the master path,
        #else return -1.
        A = line_variables_list[0]
        B = line_variables_list[1]
        C = line_variables_list[2]
        l_list = line_variables_list

        low  = 0.; high = 1.
 
        for i in range(1000):

            mid    = 0.5*(low+high)
            p_low  = self.AB.eval(low)#Node(self.AB.eval(low),"")
            p_high = self.AB.eval(high)#Node(self.AB.eval(high),"")
            p_mid  = self.AB.eval(mid)#Node(self.AB.eval(mid),"")

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
        # 此函数用来计算两点之间的距离

        d = ( (point0.x - point1.x)**2 + (point0.y - point1.y)**2 )**0.5
        return d

    def find_point(self):
        # 此函数使用二分法带返回两条任意线之间，在主线上的到次线距离为d的点。

        start = 0.0
        delta_i = 0.001
        tangent_i = 1e-7 # use a small increase ratio to get the second point to determin the tangent
        end = start + delta_i

        for num in range(1000):

            mid = 0.5*(start + end )

            C0  = self.CD.eval(start)#Node(self.CD.eval(start))
            #print "看这-->", C0, self.CD.eval(start),self.CD.eval(start).x
            
            C0_temp = self.CD.eval(start + tangent_i)

            k_C0_v   = self.get_v_tangent(C0,C0_temp)
            C0_master_r = self.lines_intersection(self.p_k_line(C0,k_C0_v))

            C1  = self.CD.eval(end)#Node(self.CD.eval(end))
            C1_temp =  self.CD.eval(end + tangent_i)# Node(self.CD.eval(end + tangent_i))
            k_C1_v   = self.get_v_tangent(C1,C1_temp)           
            C1_master_r = self.lines_intersection(self.p_k_line(C1,k_C1_v))


            Cm  = self.CD.eval(mid)#Node(self.CD.eval(mid))
            Cm_temp = self.CD.eval(mid + tangent_i)#Node(self.CD.eval(mid + tangent_i))
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

                    Ct = self.CD.eval(mid_temp)#Node(self.CD.eval(mid_temp))
                    Ct_temp = self.CD.eval(mid_temp + tangent_i)#Node(self.CD.eval(mid_temp + tangent_i))
                    k_Ct_v = self.get_v_tangent(Ct,Ct_temp)
                    Ct_master_r = self.lines_intersection(self.p_k_line(Ct,k_Ct_v))

                    if Ct_master_r == -1:
                        start_temp = mid_temp


                    elif Ct_master_r > self.tolerace:
                        end_temp = mid_temp
                        #print "start_temp2,end_temp,Ct_master_r",start_temp,end_temp,Ct_master_r
                    else:
                        start_temp = mid_temp
                        break
                    j = j+1
                    #print "j ====",j
                    if j == 500:
                        print "The loop in find_point does not reach to converge",1/0

                start = start_temp 
                # start then become the ratio with whom the vertical line could pass the most begining point
                # on master path
                end = end + delta_i

            else:

                M0 = self.AB.eval(C0_master_r)#Node(self.AB.eval(C0_master_r))
                #k_C0 = -1/k_C0_v
                #C0_M0=self.p_l_distance(M0,self.p_k_line(C0,k_C0))
                C0_M0 = self.two_p_distance(M0,C0)

                Mm = self.AB.eval(Cm_master_r)#Node(self.AB.eval(Cm_master_r))
                #k_Cm = -1/k_Cm_v
                #Cm_Mm = self.p_l_distance(Mm,self.p_k_line(Cm,k_Cm))
                Cm_Mm = self.two_p_distance(Mm,Cm)
 
                M1 = self.AB.eval(C1_master_r)#Node(self.AB.eval(C1_master_r))
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

                        Ct = self.CD.eval(mid_temp)#Node(self.CD.eval(mid_temp))
                        Ct_temp = self.CD.eval(mid_temp + tangent_i)#Node(self.CD.eval(mid_temp + tangent_i))
                        k_Ct_v = self.get_v_tangent(Ct,Ct_temp)
                        Ct_master_r = self.lines_intersection(self.p_k_line(Ct,k_Ct_v))
                        Mt = self.AB.eval(Ct_master_r)#Node(self.AB.eval(Ct_master_r))
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


