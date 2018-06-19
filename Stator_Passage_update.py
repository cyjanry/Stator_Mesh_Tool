#!/usr/bin/env python
## \Stator_Passage.py
#
""" 
Script to create a structured mesh for an Radial inflow turbine stator passage. 
The geometry and grid is defined in this file

Author: Ingo Jahn
Last modified: 22/10/2015 
"""

import numpy as np
from Stator_Profile_update import * 


   
####################################
### Setting up Basic Information ###
####################################
# For grid development, set gdata.dimensions = 2, this will create the 2-D projection of the mesh.
gdata.dimensions = 2
gdata.axisymmetric_flag = 0

# Properties are set in Coupled_Rotor.py
initial = initial_stator 
stagnation = stagnation_stator
 

# creating a trans-sonic mesh
Mesh = supersonic_mesh(R_out,R_throat,R_in,R_exit,R1,R2,throat,alpha_trailing,N_STATOR,throat_out,bl_ratio,C_ratio,corner_flag=corner_flag,r_chamfer=r_chamfer,throat_out_ratio=throat_out_ratio,theta_0=theta_0,Z = Height0 )
Mesh.calc_tuning_nodes1(Xs_ratio,I0_ratio,I1_ratio,ITOP_ratio = 0.65 )
Mesh.Line_XR1t_YR2t(X1s_ratio,TOP_ratio)
Mesh.Line_YCt2_Xt1(theta0,theta1,L1=0.1,L2=0.05,L3=D_position,Node_D="None",Node_E="None")
Mesh.Line_YCt1_Yt2(theta2,theta3,L1=0.2,L2=0.2,Node_D="None",Node_E="None")
Nose = Mesh.nose_path
Nose_BL = Mesh.nose_path_bl
Rotated_Nose = Mesh.rotated_nose_path
Rotated_Nose_BL = Mesh.rotated_nose_path_bl

Mesh.calc_depend_nodes()

# Create Lines & Block
Center = Node(0.,0.,Height0,label="Center")


centered_type = "channel"

exit_block    = "no"# "yes" # "no"
Nout = 20

# L0_BL
L0_BL_W = Line(Mesh.Xm,Mesh.Xm_BL)
L0_BL_S = Polyline([Mesh.Line_YCt2_Xt1],"",(1.-Mesh.l_ratio_Xm),1.) 
L0_BL_E = Line(Mesh.Xt1,Mesh.Xt1_BL)
L0_BL_N = Polyline([Mesh.Line_YCt2BL_Xt1BL],"",(1.-Mesh.l_ratio_Xm),1.)

# L0
L0_W = Line(Mesh.Xm_BL,Mesh.XB)
L0_N = Line(Mesh.XB,Mesh.Xt)
L0_S = L0_BL_N
L0_E = Line(Mesh.Xt1_BL,Mesh.Xt)

# L1_BL
L1_BL_W = Line(Mesh.Xt1,Mesh.Xt1_BL)
L1_BL_S = Polyline([Nose],"",0.,X1s_ratio)
L1_BL_N = Mesh.Xt1BL_XR1sBL
L1_BL_E = Line(Mesh.XR1s,Mesh.XR1s_BL)

# L1
L1_N = Line(Mesh.Xt,Mesh.Xs)
L1_W = Line(Mesh.Xt1_BL,Mesh.Xt)
L1_S = L1_BL_N
L1_E = Line(Mesh.XR1s_BL,Mesh.Xs)

# L2
L2_N = Line(Mesh.Xs,Mesh.Xi)
L2_W = Line(Mesh.XR1s_BL,Mesh.Xs)
L2_S = Line(Mesh.XR1s_BL,Mesh.I0)#Line(L1_S.eval(1.),Mesh.I0)
L2_E = Arc(Mesh.I0,Mesh.Xi,Center)


# TOPL_BL
TOPL_BL_N = Line(Mesh.XR1s,Mesh.XR1s_BL)
TOPL_BL_E = Mesh.TOPBL_XR1sBL
TOPL_BL_S = Line(Mesh.TOP,Mesh.TOP_BL)
TOPL_BL_W = Polyline([Nose],"",TOP_ratio,X1s_ratio)


# TOPL
TOPL_W = TOPL_BL_E
TOPL_E = Arc(Mesh.ITOP,Mesh.I0,Center)
TOPL_N = Line(TOPL_W.eval(1.),Mesh.I0)
TOPL_S = Line(TOPL_W.eval(0.),Mesh.ITOP) 


# TOPR_BL
if centered_type == "blade":
  TOPR_BL_N = TOPL_BL_S
  TOPR_BL_E = Polyline([Nose_BL],"",Y2s_ratio,TOP_ratio)
  TOPR_BL_W = Polyline([Nose],"",Y2s_ratio,TOP_ratio)
  TOPR_BL_S = Line(TOPR_BL_W.eval(0.),TOPR_BL_E.eval(0.))

elif centered_type == "channel": 
  TOPR_BL_N = Line(Mesh.WTOP,Mesh.WTOP_BL)# Line(Mesh.WTOP,Mesh.WTOP_BL)
  TOPR_BL_E = Polyline([Rotated_Nose_BL],"",Y2s_ratio,TOP_ratio)
  TOPR_BL_W = Polyline([Rotated_Nose],"",Y2s_ratio,TOP_ratio)
  TOPR_BL_S = Line(TOPR_BL_W.eval(0.),TOPR_BL_E.eval(0.))
  print "--",Mesh.WTOP_BL, Mesh.rotated_nose_path_bl.eval(TOP_ratio)



# TOPR 
if centered_type == "blade":
  TOPR_W = TOPR_BL_E
  TOPR_E = Arc(Mesh.I1,Mesh.ITOP,Center)
  TOPR_N = Line(TOPR_W.eval(1.),Mesh.ITOP)
  TOPR_S = Line(TOPR_W.eval(0.),Mesh.I1)  

elif centered_type == "channel": 
  
  TOPR_W = TOPR_BL_E
  TOPR_E = Arc(Mesh.I1.copy().rotate_about_zaxis(2.*np.pi/N_STATOR),Mesh.ITOP.copy().rotate_about_zaxis(2.*np.pi/N_STATOR),Center)
  TOPR_N = Line(TOPR_W.eval(1.),Mesh.ITOP.copy().rotate_about_zaxis(2.*np.pi/N_STATOR))
  TOPR_S = Line(TOPR_W.eval(0.),Mesh.I1.copy().rotate_about_zaxis(2.*np.pi/N_STATOR))


# U0_BL 
if centered_type == "blade":
  U0_BL_W = Line(Mesh.YCt1_BL,Mesh.YCt1)
  U0_BL_N = Mesh.Line_YCt1_Yt2
  U0_BL_E = Line(Mesh.Yt2_BL,Mesh.Yt2) 
  U0_BL_S = Mesh.Line_YCt1BL_Yt2BL
elif centered_type == "channel": 
  U0_BL_W = Line(Mesh.XCt1_BL,Mesh.XCt1)
  U0_BL_N = Mesh.Line_XCt1_Xt2
  U0_BL_E = Line(Mesh.Xt2_BL,Mesh.Xt2) 
  U0_BL_S = Mesh.Line_XCt1BL_Xt2BL



# U0
if centered_type == "blade":
  U0_W = Line(Mesh.YB, Mesh.YCt1_BL)
  U0_N = U0_BL_S
  U0_E = Line(Mesh.Yt, Mesh.Yt2_BL)
  U0_S = Line(Mesh.YB, Mesh.Yt)
elif centered_type == "channel": 
  U0_W = Line(Mesh.XB, Mesh.XCt1_BL)
  U0_N = U0_BL_S
  U0_E = Line(Mesh.Xt, Mesh.Xt2_BL)
  U0_S = Line(Mesh.XB, Mesh.Xt)



# U1_BL
if centered_type == "blade":
  U1_BL_W = Line(Mesh.Yt2_BL,Mesh.Yt2)#Mesh.Yt2bl_Yt2
  U1_BL_N = Polyline([Nose],"",1.,Y2s_ratio)
  U1_BL_S = Polyline([Nose_BL],"",1.,Y2s_ratio)
  U1_BL_E = Line(U1_BL_S.eval(1.), U1_BL_N.eval(1.))
elif centered_type == "channel": 
  U1_BL_W = Line(Mesh.Xt2_BL,Mesh.Xt2)#Mesh.Yt2bl_Yt2
  U1_BL_N = Polyline([Rotated_Nose],"",1.,Y2s_ratio)
  U1_BL_S = Polyline([Rotated_Nose_BL],"",1.,Y2s_ratio)
  U1_BL_E = Line(U1_BL_S.eval(1.), U1_BL_N.eval(1.))

# U1
if centered_type == "blade":
  U1_S = Line(Mesh.Yt,Mesh.Ys)
  U1_W = Line(Mesh.Yt,Mesh.Yt2_BL)
  U1_N = U1_BL_S
  U1_E = Line(Mesh.Ys,U1_N.eval(1.))
elif centered_type == "channel": 
  U1_S = Line(Mesh.Xt,Mesh.Xs)
  U1_W = Line(Mesh.Xt,Mesh.Xt2_BL)
  U1_N = U1_BL_S
  U1_E = Line(Mesh.Xs,U1_N.eval(1.))



# U2
if centered_type == "blade":
  U2_S = Line(Mesh.Ys,Mesh.Yi)
  U2_W = U1_E
  U2_N = Line(U1_N.eval(1.),Mesh.I1)
  U2_E = Arc(Mesh.Yi,Mesh.I1,Center)
elif centered_type == "channel":
  U2_S = Line(Mesh.Xs,Mesh.Xi)
  U2_W = U1_E
  U2_N = Line(U1_N.eval(1.),Mesh.I1.copy().rotate_about_zaxis(2.*np.pi/N_STATOR))
  U2_E = Arc(Mesh.Xi,Mesh.I1.copy().rotate_about_zaxis(2.*np.pi/N_STATOR),Center)


LE_BL_W = Line(Mesh.Xn,Mesh.Xn_BL)
LE_BL_N = Polyline([Mesh.Line_YCt2BL_Xt1BL],"",(1.-Mesh.l_ratio_Xn),(1.-Mesh.l_ratio_Xa1),1)
LE_BL_E = Line(Mesh.Xa1,Mesh.Xa1_BL)
LE_BL_S = Polyline([Mesh.Line_YCt2_Xt1],"",(1-Mesh.l_ratio_Xn),(1.-Mesh.l_ratio_Xa1),1)


#LT0_BL_W = Line(Mesh.YCt0_BL,Mesh.YCt0)
#LT0_BL_N = Arc(Mesh.YCt0,Mesh.YCt1,Mesh.YC)
#LT0_BL_E = U0_BL_W
#LT0_BL_S = Arc(Mesh.YCt0_BL,Mesh.YCt1_BL,Mesh.YC)


LT2_BL_W = Line(Mesh.YCt2,Mesh.YCt2_BL)
LT2_BL_N = Polyline([Mesh.Line_YCt2BL_Xt1BL],"",0.,(1.-Mesh.l_ratio_Xn),1)
LT2_BL_E = Line(Mesh.Xn,Mesh.Xn_BL)
LT2_BL_S = Polyline([Mesh.Line_YCt2_Xt1],"",0.,(1-Mesh.l_ratio_Xn),1)



if centered_type == "blade":
  LT1_BL_W = Line(Mesh.YCt0,Mesh.YCt0_BL)
  LT1_BL_N = Arc(Mesh.YCt0_BL,Mesh.YCt2_BL,Mesh.YC)
  LT1_BL_E = LT2_BL_W
  LT1_BL_S = Arc(Mesh.YCt0,Mesh.YCt2,Mesh.YC)
elif centered_type == "channel":
  LT1_BL_W = Line(Mesh.XCt0,Mesh.XCt0_BL)
  LT1_BL_N = Arc(Mesh.XCt0_BL,Mesh.XCt2_BL,Mesh.XC)
  LT1_BL_E = Line(Mesh.XCt2,Mesh.XCt2_BL)
  LT1_BL_S = Arc(Mesh.XCt0,Mesh.XCt2,Mesh.XC)


if centered_type == "blade":
  LT1_W = Line(Mesh.YCt0_BL,Mesh.YCt0e)#Bezier([Mesh.YCt0_BL,Mesh.YCt0_ctrl,Mesh.YCt0e])#Mesh.Line_YCt0BL_YCt0e
  LT1_N = Arc(Mesh.YCt0e,Mesh.YCte,Center)
  LT1_E = Line(Mesh.YCt2_BL,Mesh.YCte)
  LT1_S = Arc(Mesh.YCt0_BL,Mesh.YCt2_BL,Mesh.YC)
elif centered_type == "channel":
  LT1_W = Line(Mesh.XCt0_BL,Mesh.XCt0e)#Bezier([Mesh.YCt0_BL,Mesh.YCt0_ctrl,Mesh.YCt0e])#Mesh.Line_YCt0BL_YCt0e
  LT1_N = Arc(Mesh.XCt0e,Mesh.XCte,Center)
  LT1_E = Line(Mesh.XCt2_BL,Mesh.XCte)
  LT1_S = Arc(Mesh.XCt0_BL,Mesh.XCt2_BL,Mesh.XC)  

LT2_W = Line(Mesh.YCt2_BL,Mesh.YCte)
LT2_N = Arc(Mesh.YCte,Mesh.XeL,Center)
LT2_E = Line(Mesh.Xn_BL,Mesh.XeL)
LT2_S = LT2_BL_N

LE_W = Line(Mesh.Xn_BL,Mesh.XeL)
LE_N = Polyline([Mesh.Line_XeL_X,Line(Mesh.X,Mesh.Xa)],"",0,1.,1)#  Line(Mesh.Xe,Mesh.Xa)
LE_E = Line(Mesh.Xa1_BL,Mesh.Xa)
LE_S = Polyline([Mesh.Line_YCt2BL_Xt1BL],"",(1.-Mesh.l_ratio_Xn),(1-Mesh.l_ratio_Xa1),1)



LM_BL_W = Line(Mesh.Xa1,Mesh.Xa1_BL)
LM_BL_N = Polyline([Mesh.Line_YCt2BL_Xt1BL],"",(1.-Mesh.l_ratio_Xa1),(1.-Mesh.l_ratio_Xm),1)
LM_BL_E = Line(Mesh.Xm,Mesh.Xm_BL)
LM_BL_S = Polyline([Mesh.Line_YCt2_Xt1],"",(1.-Mesh.l_ratio_Xa1),(1-Mesh.l_ratio_Xm),1)
 
LM_W    = Line(Mesh.Xa1_BL, Mesh.Xa)
LM_N    = Line(Mesh.Xa, Mesh.XB)
LM_E    = Line(Mesh.Xm_BL, Mesh.XB)
LM_S    = LM_BL_N



if centered_type == "blade":
  UM_BL_W = Line(Mesh.YCt0_BL,Mesh.YCt0)
  UM_BL_N = Arc(Mesh.YCt0,Mesh.YCt1,Mesh.YC)
  UM_BL_E = Line(Mesh.YCt1_BL,Mesh.YCt1)
  UM_BL_S = Arc(Mesh.YCt0_BL, Mesh.YCt1_BL,Mesh.YC)
elif centered_type == "channel":
  UM_BL_W = Line(Mesh.XCt0_BL,Mesh.XCt0)
  UM_BL_N = Arc(Mesh.XCt0,Mesh.XCt1,Mesh.XC)
  UM_BL_E = Line(Mesh.XCt1_BL,Mesh.XCt1)
  UM_BL_S = Arc(Mesh.XCt0_BL, Mesh.XCt1_BL,Mesh.XC)

if centered_type == "blade":
  UM_W = Bezier([Mesh.Ya,Mesh.YCt0_ctrl,Mesh.YCt0_BL],"",0.,1.,1)
  UM_N = UM_BL_S
  UM_E = U0_W
  UM_S = Line(Mesh.Ya,Mesh.YB)
elif centered_type == "channel":
  UM_W = Bezier([Mesh.Xa,Mesh.XCt0_ctrl,Mesh.XCt0_BL],"",0.,1.,1)
  UM_N = UM_BL_S
  UM_E = U0_W
  UM_S = Line(Mesh.Xa,Mesh.XB)



if centered_type == "blade":
  UT_W = Arc(Mesh.YeL,Mesh.YCt0e,Center)
  UT_N = Polyline([LT1_W],"",1.0,0.0,1)
  UT_E = UM_W
  UT_S = Polyline([Mesh.Line_YeL_Y,Line(Mesh.Y,Mesh.Ya)],"",0.,1.,1)
elif centered_type == "channel":
  UT_W = Arc(Mesh.XeL,Mesh.XCt0e,Center)
  UT_N = Polyline([Line(Mesh.XCt0_BL,Mesh.XCt0e)],"",1.0,0.0,1)
  UT_E = UM_W
  UT_S = Polyline([Mesh.Line_XeL_X,Line(Mesh.X,Mesh.Xa)],"",0.,1.,1)



if exit_block == "yes" : 
  if centered_type == "channel":
    LT1O_W = Arc(Mesh.XCt0e_O,Mesh.XCte_O,Center)
    LT1O_N = Line(Mesh.XCte_O,Mesh.XCte)
    LT1O_E = Arc(Mesh.XCt0e,Mesh.XCte,Center)
    LT1O_S = Line(Mesh.XCt0e_O,Mesh.XCt0e)

    UTO_W  = Arc(Mesh.XeL_O, Mesh.XCt0e_O, Center)
    UTO_N  = Line(Mesh.XCt0e_O,Mesh.XCt0e)
    UTO_E  = Arc(Mesh.XeL, Mesh.XCt0e, Center)
    UTO_S  = Line(Mesh.XeL_O,Mesh.XeL)

    LT2O_W = Arc(Mesh.YCte_O,Mesh.XeL_O,Center)
    LT2O_N = Line(Mesh.XeL_O,Mesh.XeL)
    LT2O_E = Arc(Mesh.YCte,Mesh.XeL,Center)
    LT2O_S = Line(Mesh.YCte_O,Mesh.YCte)    








#0-4


print "creating L2 ..."
L2 = Block2D(make_patch(L2_N, L2_E, L2_S, L2_W), nni=Nin, nnj=Nlow,
     cf_list=[S_t2,None,L2_TOP,L2_L1],fill_condition=initial, label="L2")



print "creating L1_BL ..."
L1_BL = Block2D(make_patch(L1_BL_N, L1_BL_E, L1_BL_S, L1_BL_W), nni=Nnoz1, nnj=Nbl,
      cf_list=[No_3,BL_l,No_4,BL_l],fill_condition=initial, label="L1_BL")

print "creating L1 ..."
L1 = Block2D(make_patch(L1_N, L1_E, L1_S, L1_W), nni=Nnoz1, nnj=Nlow,
     cf_list=[No_2,L2_L1,No_3,None],fill_condition=initial, label="L1")

print "creating L0_BL ..."
L0_BL = Block2D(make_patch(L0_BL_N, L0_BL_E, L0_BL_S, L0_BL_W), nni=Nnoz0, nnj=Nbl,
     cf_list=[No_30,BL_l,No_40,BL_l], fill_condition=initial, label="L0_BL")
print "creating L0 ..."
L0 = Block2D(make_patch(L0_N, L0_E, L0_S, L0_W), nni=Nnoz0, nnj=Nlow,
     cf_list=[No_20,None,No_30,None], fill_condition=initial, label="L0")



#5-9
print "creating TOPL ..."
#if   centered_type == "channel":
#    data = make_patch(TOPL_N, TOPL_E, TOPL_S, TOPL_W).clone().rotate_about_zaxis(2.*np.pi/N_STATOR)
#elif centered_type == "blade":
data = make_patch(TOPL_N, TOPL_E, TOPL_S, TOPL_W)
TOPL = Block2D(data, nni=Nin, nnj=Ntopl,
     cf_list=[L2_TOP,None,TOPR_TOPL,S_t1],fill_condition=initial, label="TOP0")



print "creating TOPL_BL ..."
TOPL_BL = Block2D(make_patch(TOPL_BL_N, TOPL_BL_E, TOPL_BL_S, TOPL_BL_W), nni=Nbl, nnj=Ntopl,
     cf_list=[BL_l,S_t1,BL_l,S_t1],fill_condition=initial, label="TOPL_BL")






print "creating TOPR ..."
data = make_patch(TOPR_N, TOPR_E, TOPR_S, TOPR_W)

TOPR = Block2D(data, nni=Nin, nnj=Ntopr,
     cf_list=[TOPR_TOPL,None,S_t0,S_t1],fill_condition=initial, label="TOP0")


print "creating TOPR_BL ..."
data = make_patch(TOPR_BL_N, TOPR_BL_E, TOPR_BL_S, TOPR_BL_W)
TOPR_BL = Block2D(data, nni=Nbl, nnj=Ntopr,
     cf_list=[BL_l,S_t1,BL_l,S_t1],fill_condition=initial, label="TOPR_BL")




print "creating U2 ..."

data = make_patch(U2_N, U2_E, U2_S, U2_W)

U2 = Block2D(data, nni=Nin, nnj=Nup,
     cf_list=[S_t0,None,S_t2,None],fill_condition=initial, label="U2")



#10-14
print "creating U1_BL ..."
#if   centered_type == "channel":
#  data = make_patch(U1_BL_N, U1_BL_E, U1_BL_S, U1_BL_W).clone().rotate_about_zaxis(2.*np.pi/N_STATOR)
#elif  centered_type == "blade":
data =  make_patch(U1_BL_N, U1_BL_E, U1_BL_S, U1_BL_W)
U1_BL = Block2D(data, nni=Nnoz1, nnj=Nbl,
     cf_list=[No_0,BL_r,No_1,BL_r],fill_condition=initial, label="U1_BL")



print "creating U1 ..."
data = make_patch(U1_N, U1_E, U1_S, U1_W)
U1 = Block2D(data, nni=Nnoz1, nnj=Nup,
     cf_list=[No_1,None,No_2,None],fill_condition=initial, label="U1")


print "creating U0_BL ..."
#if   centered_type == "channel":
#  data = make_patch(U0_BL_N, U0_BL_E, U0_BL_S, U0_BL_W).clone().rotate_about_zaxis(2.*np.pi/N_STATOR)
#elif centered_type == "blade":
data = make_patch(U0_BL_N, U0_BL_E, U0_BL_S, U0_BL_W)
U0_BL = Block2D(data, nni=Nnoz0, nnj=Nbl,
     cf_list=[No_00,BL_r,No_10,BL_c],fill_condition=initial, label="U0_BL")

print "creating U0 ..."
#if centered_type == "channel":
#  data = make_patch(U0_N, U0_E, U0_S, U0_W).clone().rotate_about_zaxis(2.*np.pi/N_STATOR)
#elif  centered_type == "blade":
data = make_patch(U0_N, U0_E, U0_S, U0_W)
U0 = Block2D(data, nni=Nnoz0, nnj=Nup,
     cf_list=[No_10,None,No_20,UM_U0],fill_condition=initial, label="U0")



print "creating LM_BL ..."
LM_BL = Block2D(make_patch(LM_BL_N, LM_BL_E, LM_BL_S, LM_BL_W), nni=Nnozm, nnj=Nbl,
     cf_list=[None,BL_l,None,BL_l],fill_condition=initial, label="LM_BL")




print "creating LM ..."
LM = Block2D(make_patch(LM_N, LM_E, LM_S, LM_W), nni=Nnozm, nnj=Nlow,
     cf_list=[None,None,None,None],fill_condition=initial, label="LM")


#15-19
print "creating LE_BL ..."
LE_BL = Block2D(make_patch(LE_BL_N, LE_BL_E, LE_BL_S, LE_BL_W), nni=Nlow, nnj=Nbl,
     cf_list=[LE_LEBL,BL_l,LE_LEBL,BL_l],fill_condition=initial, label="LE_BL")



print "creating LE ..."
LE = Block2D(make_patch(LE_N, LE_E, LE_S, LE_W), nni=Nlow, nnj=Nlow,
     cf_list=[LE_UT,None,LE_LEBL,LE_LT2],fill_condition=initial, label="LE")







print "creating LT2_BL ..."
LT2_BL = Block2D(make_patch(LT2_BL_N, LT2_BL_E, LT2_BL_S, LT2_BL_W), nni=Nnozt2, nnj=Nbl,
     cf_list=[None,BL_l,None,BL_l],fill_condition=initial, label="LT2_BL")

#20-24
print "creating LT2 ..."
LT2 = Block2D(make_patch(LT2_N, LT2_E, LT2_S, LT2_W), nni=Nnozt2, nnj=Nlow,
     cf_list=[None,LE_LT2,None,LT1_LT2],fill_condition=initial, label="LT2")

print "creating LT1_BL ..."
LT1_BL = Block2D(make_patch(LT1_BL_N, LT1_BL_E, LT1_BL_S, LT1_BL_W), nni=Nnozt1, nnj=Nbl,
     cf_list=[None,BL_l,None,BL_l],fill_condition=initial, label="LT1_BL")

#print "creating LT0_BL ..."
#LT0_BL = Block2D(make_patch(LT0_BL_N, LT0_BL_E, LT0_BL_S, LT0_BL_W), nni=Nnoz0, nnj=Nupp,
#     cf_list=[No_10,None,No_20,M1_U0],fill_condition=initial, label="LT0_BL")


print "creating LT1 ..."
LT1 = Block2D(make_patch(LT1_N, LT1_E, LT1_S, LT1_W), nni=Nnozt1, nnj=Nlow,
     cf_list=[None,LT1_LT2,None,LT1_UT],fill_condition=initial, label="LT1")

print "creating UM_BL ..."
data = make_patch(UM_BL_N, UM_BL_E, UM_BL_S, UM_BL_W)
UM_BL = Block2D(data, nni=Nnozm, nnj=Nbl,
     cf_list=[None,BL_r,None,BL_r],fill_condition=initial, label="UM_BL")

print "creating UM ..."
data = make_patch(UM_N, UM_E, UM_S, UM_W)
UM = Block2D(data, nni=Nnozm, nnj=Nup,
     cf_list=[None,UM_U0,None,None],fill_condition=initial, label="UM")


print "creating UT ..."
data = make_patch(UT_N, UT_E, UT_S, UT_W)
UT = Block2D(data, nni=Nlow, nnj=Nup,
     cf_list=[UT_LT1,None,LE_UT,None],fill_condition=initial, label="UT")


if exit_block == "yes":
  print "creating LT1O ..."
  data = make_patch(LT1O_N, LT1O_E, LT1O_S, LT1O_W)
  LT1O = Block2D(data, nni=Nout, nnj=Nnozt1,
     cf_list=[UT_LT1,None,LE_UT,None],fill_condition=initial, label="LT1O")

  print "creating UTO  ..."
  data = make_patch(UTO_N, UTO_E, UTO_S, UTO_W)
  UTO = Block2D(data, nni=Nout, nnj=Nup,
     cf_list=[UT_LT1,None,LE_UT,None],fill_condition=initial, label="UTO")

  print "creating LT2O ..."
  data = make_patch(LT2O_N, LT2O_E, LT2O_S, LT2O_W)
  LT2O = Block2D(data, nni=Nout, nnj=Nnozt2,
     cf_list=[UT_LT1,None,LE_UT,None],fill_condition=initial, label="LT2O")
# link blocks
identify_block_connections()




####################################
###  define B/C                  ###
####################################
"""
OF_inlet_00 --> Inlet
OF_outlet_00 --> Outlet

OF_wall_00 --> Stator Vane
OF_wall_01 --> Vane Trailing Edge
OF_wall_02 --> Top Wall
OF_wall_03 --> Bottom Wall
"""
# set correct labels for Periodic Boundaries.
if centered_type == "blade":
  L2.bc_list[NORTH]    = SlipWallBC(label='OF_symmetry_00')
  L1.bc_list[NORTH]    = SlipWallBC(label='OF_symmetry_00')
  L0.bc_list[NORTH]    = SlipWallBC(label='OF_symmetry_00')
  LM.bc_list[NORTH]    = SlipWallBC(label='OF_symmetry_00')
  LE.bc_list[NORTH]    = SlipWallBC(label='OF_symmetry_00')


  U2.bc_list[SOUTH]    = SlipWallBC(label='OF_symmetry_01')
  U1.bc_list[SOUTH]    = SlipWallBC(label='OF_symmetry_01')
  U0.bc_list[SOUTH]    = SlipWallBC(label='OF_symmetry_01')
  UM.bc_list[SOUTH]    =  SlipWallBC(label='OF_symmetry_01') 
  UT.bc_list[SOUTH]    =  SlipWallBC(label='OF_symmetry_01')




elif centered_type == "channel":

  TOPL.bc_list[SOUTH]       = SlipWallBC(label='OF_symmetry_00')
  TOPL_BL.bc_list[SOUTH]    = SlipWallBC(label='OF_symmetry_00')
  LT2_BL.bc_list[WEST]      = SlipWallBC(label='OF_symmetry_00')
  LT2.bc_list[WEST]         = SlipWallBC(label='OF_symmetry_00')

  TOPR.bc_list[NORTH]       = SlipWallBC(label='OF_symmetry_01')
  TOPR_BL.bc_list[NORTH]    = SlipWallBC(label='OF_symmetry_01')
  LT1_BL.bc_list[EAST]      = SlipWallBC(label='OF_symmetry_01') 
  LT1.bc_list[EAST]         = SlipWallBC(label='OF_symmetry_01')



L2.bc_list[EAST] = SubsonicInBC(stagnation,label='OF_inlet_00',
                    direction_type='radial',direction_alpha=inlet_flow_direction)
TOPL.bc_list[EAST] = SubsonicInBC(stagnation,label='OF_inlet_00',
                    direction_type='radial',direction_alpha=inlet_flow_direction)
TOPR.bc_list[EAST] = SubsonicInBC(stagnation,label='OF_inlet_00',
                    direction_type='radial',direction_alpha=inlet_flow_direction)
U2.bc_list[EAST] = SubsonicInBC(stagnation,label='OF_inlet_00',
                    direction_type='radial',direction_alpha=inlet_flow_direction)


TOPL_BL.bc_list[WEST] = AdiabaticBC(label='OF_wall_00')
L1_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
L0_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
LM_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
LE_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
LT2_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')


TOPR_BL.bc_list[WEST] = AdiabaticBC(label='OF_wall_01')
U1_BL.bc_list[NORTH] = AdiabaticBC(label='OF_wall_01')
U0_BL.bc_list[NORTH] = AdiabaticBC(label='OF_wall_01')
UM_BL.bc_list[NORTH] = AdiabaticBC(label='OF_wall_01')
LT1_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_01')



if exit_block == "yes":
  LT1O.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_01')
  LT2O.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_00')
  
  LT1O.bc_list[WEST]  = FixedPOutBC(Pout, label='OF_outlet_00')
  UTO.bc_list[WEST]   =  FixedPOutBC(Pout, label='OF_outlet_00')
  LT2O.bc_list[WEST]  = FixedPOutBC(Pout, label='OF_outlet_00')


else:

  LT1.bc_list[NORTH] = FixedPOutBC(Pout, label='OF_outlet_00')
  LT2.bc_list[NORTH] = FixedPOutBC(Pout, label='OF_outlet_00')
  UT.bc_list[WEST]   = FixedPOutBC(Pout, label='OF_outlet_00') 
  
    
#Coords = Mesh.Xi
#r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
#T1_in = t
#Coords = Mesh.Yi
#r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
#T0_in = t


#Coords = NGV_out.E0_EX.eval(0.,1.)
#r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
#T1_out = t
#Coords = NGV_out.E3_EX.eval(0.,0.)
#r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
#T0_out = t

'''
print "####################################################"
print "##           Settings for Rotor meshing           ##"
print "####################################################"
print ""
print "Stator Inlet:" 
print "Radius            : ", R_in, " [m]"
print "Arc start     : ", T1_in, " [rad]"
print "Arc end       : ", T0_in, " [rad]"
print "Arc           : ", T1_in-T0_in, " [rad]"
print ""
print "Stator Outlet:" 
print "Radius            : ", R_exit, " [m]"
print "Z-Position Hub    : ", Height0, " [m]"
print "Z-Position Shroud : ", Height1, " [m]"  
print "Arc start         : ", T1_out, " [rad]"
print "Arc end           : ", T0_out, " [rad]"
print "Arc               : ", T1_out-T0_out, " [rad]"
print ""
print ""
print "Number of Stators : ", N_STATOR
print "####################################################"
'''

sketch.prefer_bc_labels_on_faces() # required to allow grouping of boundaries by e3prepToFoam.py

if gdata.dimensions == 2:
    # This is to make a nice *.svg file of the 2-D projection of the mesh
   #sketch.xaxis(0.05, 0.05, 0.02, 0.0)
   #sketch.yaxis(0.05, 0.05, 0.02, 0.0)
   # sketch.window(0., -0.025, 0.05, 0.025, 0.01, 0.01, 0.2, 0.2) 
   #sketch.window(xmin=0.065, ymin=0.002, xmax=0.064, ymax=0.0025)

   sketch.window(xmin=0.01, ymin=-0.005, xmax=0.045, ymax=0.030, )   #focus on the whole stator for drawing
   #sketch.window(xmin=0.022, ymin=0.005, xmax=0.023, ymax=0.006) #focus on trailing edge

   #sketch.window(xmin=0.063, ymin=-0.0015, xmax=0.066, ymax=0.0015) # Focus on trailing edge

   #sketch.do_labels(label_state=False)  # use this script to silence the labels






