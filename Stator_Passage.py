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
from Stator_Profile import * 

print R_out, corner_flag,d_chamfer




if mesh_type is "transsonic":    
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
    Mesh = supersonic_mesh(R_out,R_throat,R_in,R1,R2,throat,alpha_trailing,N_STATOR,throat_out,bl_ratio,corner_flag=corner_flag,d_chamfer=d_chamfer,throat_out_ratio=throat_out_ratio,theta_0=theta_0,Z = Height0)
    Mesh.calc_fixed_nodes1()
    Mesh.calc_flex_nodes1() 
    Mesh.calc_tuning_nodes1(Xs_ratio,I0_ratio,I1_ratio,ITOP_ratio = 0.65 )
    Mesh.Line_XR1t_YR2t(X1s_ratio,TOP_ratio )
    Nose = Mesh.nose_path
    Nose_BL = Mesh.nose_path_bl
    
    
    if corner_flag == 1:
        #####
        u0_bl_n_origin = Mesh.Line_YO2_Yt1(0.05,alpha2="None",L1=0.1,L2=0.1,Node_D="None",Node_E="None") # JQ 19/11/2016
        #Mesh.Line_YOr0_YOr1()      # JQ 19/11/2016
        u0_bl_n = Mesh.Line_Yt2_YOr1
        l0_bl_s = Mesh.Line_XO1_Xt1(0.05,alpha1="None",L1=0.1,L2=0.1,Node_D="None",Node_E="None")
        l0_bl_s_ratio = get_ratio(l0_bl_s,Mesh.XO3)

        line_Xm_XO3 = Mesh.Line_Xm_XO3()
        #line_Xm_XO3_raio = get_ratio(line_Xm_XO3,Mesh.XO3_BL) # get the ratio for point XO3_BL on line Line_Xm_XO3

        l0_bl_n = move_and_stretch( Mesh.XO1_BL , Mesh.Xt1_BL ,l0_bl_s).Line()
        l0_bl_n_ratio = get_ratio(l0_bl_n,Mesh.XO3_BL)
        
        u0_bl_s = Mesh.Line_YO2BL_Yt2BL#move_and_stretch( Mesh.YO2_BL , Mesh.Yt2_BL , u0_bl_n_origin).Line()
        u0_bl_s_ratio = get_ratio(u0_bl_s,Mesh.YOr1_BL)

        # Generate some Bezier curves:
        #Mesh.Find_ctrl_points_U0M1()
        Mesh.Initialize_Bezier_Curves()
        ####  JQ     

    # Create Lines & Block
    Center = Node(0.,0.,Height0,label="Center")


    if corner_flag == 1:
        # L0_BL
        L0_BL_W = Line(Mesh.XO3,Mesh.XO3_BL)#Polyline([line_Xm_XO3],"",0.,line_Xm_XO3_raio) #Arc(Mesh.XO1,Mesh.XO1_BL,Center)
        L0_BL_S = Polyline([l0_bl_s],"",l0_bl_s_ratio,1.) 
        L0_BL_E = Line(Mesh.Xt1,Mesh.Xt1_BL)
        L0_BL_N = Polyline([l0_bl_n],"",l0_bl_n_ratio,1.)#l0_bl_n #L0_BL_N_temp.Line()

        # L0
        L0_W = Mesh.L0M0#Line(Mesh.XO3_BL,Mesh.Xm)#Polyline([line_Xm_XO3],"",line_Xm_XO3_raio,1.) #Arc(Mesh.XO1_BL,Mesh.X,Center)
        L0_N = Line(Mesh.Xm,Mesh.Xt)
        L0_S = L0_BL_N
        L0_E = Line(Mesh.Xt1_BL,Mesh.Xt)

        # L1_BL
        L1_BL_W = L0_BL_E
        L1_BL_S = Polyline([Nose],"",0.,X1s_ratio)
        L1_BL_N = Mesh.Xt1BL_XR1sBL#Line(Mesh.Xt1_BL,Mesh.XR1s_BL)#Polyline([Nose_BL],"",0.,X1s_ratio)
        L1_BL_E = Line(Mesh.XR1s,Mesh.XR1s_BL)#Line(L1_BL_S.eval(1.), L1_BL_N.eval(1.))
        print Mesh.Xt1,Mesh.Xt1_BL,Mesh.XR1s_BL,Mesh.XR1s
        # L1
        L1_N = Line(Mesh.Xt,Mesh.Xs)
        L1_W = L0_E
        L1_S = L1_BL_N
        L1_E = Line(Mesh.XR1s_BL,Mesh.Xs)#Line(L1_S.eval(1.),Mesh.Xs)

        # L2
        L2_N = Line(Mesh.Xs,Mesh.Xi)
        L2_W = L1_E
        L2_S = Line(Mesh.XR1s_BL,Mesh.I0)#Line(L1_S.eval(1.),Mesh.I0)
        L2_E = Arc(Mesh.I0,Mesh.Xi,Center)

        # U0_BL 
        #U0_BL_W = Arc(Mesh.YO2_BL,Mesh.YO2,Center)
        U0_BL_W = Line(Mesh.YOr1_BL,Mesh.YOr1)#Line(Mesh.YO2_BL, Mesh.YOr1)
        U0_BL_N = u0_bl_n # changing JQ 19/11/2016
        U0_BL_E = Mesh.Yt2bl_Yt2#Line(Mesh.Yt2_BL,Mesh.Yt2) 
        #U0_BL_S_temp = move_and_stretch( Mesh.YO2_BL , Mesh.Yt2_BL , U0_BL_N)
        U0_BL_S = Polyline([u0_bl_s],"",u0_bl_s_ratio,1.)  #U0_BL_S_temp.Line()

        # U0
        U0_W = Mesh.U0M1 #Line(Mesh.Ym,Mesh.YOr1_BL) #Arc(Mesh.Y,Mesh.YO2_BL,Center)
        U0_S = Line(Mesh.Ym,Mesh.Yt)
        U0_N = U0_BL_S
        U0_E = Mesh.Yt_Yt2bl#Line(Mesh.Yt,Mesh.Yt2_BL) 
    elif corner_flag == 0:
        # L0_BL
        L0_BL_W = Arc(Mesh.XO1,Mesh.XO1_BL,Center)
        L0_BL_S = Mesh.Line_XO1_Xt1(0.05,alpha1="None",L1=0.1,L2=0.1,Node_D="None",Node_E="None")
        L0_BL_E = Line(Mesh.Xt1,Mesh.Xt1_BL)
        L0_BL_N_temp = move_and_stretch( Mesh.XO1_BL , Mesh.Xt1_BL , L0_BL_S)
        L0_BL_N = L0_BL_N_temp.Line()
     
        # L0
        L0_W = Arc(Mesh.XO1_BL,Mesh.X,Center)
        L0_N = Line(Mesh.X,Mesh.Xt)
        L0_S = L0_BL_N
        L0_E = Line(Mesh.Xt1_BL,Mesh.Xt)
 
        # L1_BL
        L1_BL_W = L0_BL_E
        L1_BL_S = Polyline([Nose],"",0.,X1s_ratio)
        L1_BL_N = Line(Mesh.Xt1_BL,Mesh.XR1s_BL)#Polyline([Nose_BL],"",0.,X1s_ratio)
        L1_BL_E = Line(Mesh.XR1s_BL,Mesh.TOP_BL)#Line(L1_BL_S.eval(1.), L1_BL_N.eval(1.))

        # L1
        L1_N = Line(Mesh.Xt,Mesh.Xs)
        L1_W = L0_E
        L1_S = L1_BL_N
        L1_E = Line(Mesh.XR1s_BL,Mesh.Xs)#Line(L1_S.eval(1.),Mesh.Xs)

        # L2
        L2_N = Line(Mesh.Xs,Mesh.Xi)
        L2_W = L1_E
        L2_S = Line(Mesh.XR1s_BL,Mesh.I0)#Line(L1_S.eval(1.),Mesh.I0)
        L2_E = Arc(Mesh.I0,Mesh.Xi,Center)

        # U0_BL 
        U0_BL_W = Arc(Mesh.YO2_BL,Mesh.YO2,Center)
        U0_BL_N = Mesh.Line_YO2_Yt1(0.05,alpha2="None",L1=0.1,L2=0.1,Node_D="None",Node_E="None")
        U0_BL_E = Line(Mesh.Yt2_BL,Mesh.Yt2)
        U0_BL_S_temp = move_and_stretch( Mesh.YO2_BL , Mesh.Yt2_BL , U0_BL_N)
        U0_BL_S = U0_BL_S_temp.Line()
     
        # U0
        U0_W = Arc(Mesh.Y,Mesh.YO2_BL,Center)
        U0_S = Line(Mesh.Y,Mesh.Yt)
        U0_N = U0_BL_S
        U0_E = Line(Mesh.Yt,Mesh.Yt2_BL)
    else:
        print "Wrong value for corner flag"

                                                                     

    # U1_BL
    U1_BL_W = U0_BL_E
    U1_BL_N = Polyline([Nose],"",1.,Y2s_ratio)
    U1_BL_S = Polyline([Nose_BL],"",1.,Y2s_ratio)
    U1_BL_E = Line(U1_BL_S.eval(1.), U1_BL_N.eval(1.))

    # U1
    U1_S = Line(Mesh.Yt,Mesh.Ys)
    U1_W = U0_E
    U1_N = U1_BL_S
    U1_E = Line(Mesh.Ys,U1_N.eval(1.))

    # U2
    U2_S = Line(Mesh.Ys,Mesh.Yi)
    U2_W = U1_E
    U2_N = Line(U1_N.eval(1.),Mesh.I1)
    U2_E = Arc(Mesh.Yi,Mesh.I1,Center)

    # TOPL_BL
    TOPL_BL_N = Line(Mesh.XR1s,Mesh.XR1s_BL)
    TOPL_BL_E = Mesh.TOPBL_XR1SBL#Line(Mesh.TOP_BL,Mesh.XR1s_BL)
    TOPL_BL_S = Line(Mesh.TOP,Mesh.TOP_BL)
    TOPL_BL_W = Polyline([Nose],"",TOP_ratio,X1s_ratio)
    #TOPL_BL_W = Polyline([Nose],"",Y2s_ratio,X1s_ratio)
    #TOPL_BL_E = Polyline([Nose_BL],"",Y2s_ratio,X1s_ratio)
    #TOPL_BL_N = Line(TOP_BL_W.eval(1.),TOP_BL_E.eval(1.))
    #TOPL_BL_S = Line(TOP_BL_W.eval(0.),TOP_BL_E.eval(0.))

    # TOPR_BL
    TOPR_BL_N = TOPL_BL_S
    TOPR_BL_E = Polyline([Nose_BL],"",Y2s_ratio,TOP_ratio)
    TOPR_BL_W = Polyline([Nose],"",Y2s_ratio,TOP_ratio)
    TOPR_BL_S = Line(TOPR_BL_W.eval(0.),TOPR_BL_E.eval(0.))
    #TOP_BL_W = Polyline([Nose],"",Y2s_ratio,X1s_ratio)
    #TOP_BL_E = Polyline([Nose_BL],"",Y2s_ratio,X1s_ratio)
    #TOP_BL_N = Line(TOP_BL_W.eval(1.),TOP_BL_E.eval(1.))
    #TOP_BL_S = Line(TOP_BL_W.eval(0.),TOP_BL_E.eval(0.))

    # TOPL
    TOPL_W = TOPL_BL_E
    TOPL_E = Arc(Mesh.ITOP,Mesh.I0,Center)
    TOPL_N = Line(TOPL_W.eval(1.),Mesh.I0)
    TOPL_S = Line(TOPL_W.eval(0.),Mesh.ITOP)   


    # TOPR 
    TOPR_W = TOPR_BL_E
    TOPR_E = Arc(Mesh.I1,Mesh.ITOP,Center)
    TOPR_N = Line(TOPR_W.eval(1.),Mesh.ITOP)
    TOPR_S = Line(TOPR_W.eval(0.),Mesh.I1)   


    # TOP_BL
#    TOP_BL_W = Polyline([Nose],"",Y2s_ratio,X1s_ratio)
#    TOP_BL_E = Polyline([Nose_BL],"",Y2s_ratio,X1s_ratio)
#    TOP_BL_N = Line(TOP_BL_W.eval(1.),TOP_BL_E.eval(1.))
#    TOP_BL_S = Line(TOP_BL_W.eval(0.),TOP_BL_E.eval(0.))

    # TOP
#    TOP_W = TOP_BL_E
#    TOP_E = Arc(Mesh.I1,Mesh.I0,Center)
#    TOP_N = Line(TOP_W.eval(1.),Mesh.I0)
#    TOP_S = Line(TOP_W.eval(0.),Mesh.I1)

    # Create the Outlet Mesh
    if corner_flag == 1:
        # Origin
        M0_L0_N = Line(Mesh.X,Mesh.Xm) 
        M1_U0_S = Line(Mesh.Y,Mesh.Ym)
        l0_w = Arc(Mesh.XO1_BL,Mesh.X,Center)
        l0_bl_w = Arc(Mesh.XO1,Mesh.XO1_BL,Center)   
        u0_bl_w = Arc(Mesh.YO2_BL,Mesh.YO2,Center)
        u0_w = Arc(Mesh.Y,Mesh.YO2_BL,Center)

        NGV_out = NGV_outlet_chamfer([Mesh.X,Mesh.XO1_BL,Mesh.XO1,Mesh.YO2,Mesh.YO2_BL,Mesh.Y,Mesh.YOr0,Mesh.YOr1,Mesh.O_chamfer,Mesh.YOr0_BL,Mesh.YOr1_BL,Mesh.XO1_BL_2,Mesh.Xm,Mesh.Ym,Mesh.XO4_BL,Mesh.XO4],[M0_L0_N,l0_bl_n,l0_bl_s,u0_bl_n_origin,u0_bl_s,M1_U0_S],[l0_w,l0_bl_w,u0_bl_w,u0_w],R_out,R_exit,bl_ratio,N_STATOR,Vertical_Line=Vertical_Line, alpha_exit=alpha_exit, arc=arc)
    elif corner_flag == 0:
        NGV_out = NGV_outlet([Mesh.X,Mesh.XO1_BL,Mesh.XO1,Mesh.YO2,Mesh.YO2_BL,Mesh.Y],
                [L0_N,L0_S,L0_BL_S,U0_BL_N,U0_N,U0_S],[L0_W,L0_BL_W,U0_BL_W,U0_W],
                R_out,R_exit, Vertical_Line=Vertical_Line, alpha_exit=alpha_exit, arc=arc)
        #NGV_out = NGV_outlet([Mesh.X,Mesh.XO1_BL,Mesh.XO1,Mesh.YO2,Mesh.YO2_BL,Mesh.Y,Mesh.YOr0,Mesh.YOr1,Mesh.O_chamfer,Mesh.YOr0_BL,Mesh.YOr1_BL,Mesh.XO1_BL_2,Mesh.Xm,Mesh.Ym,Mesh.XO4_BL,Mesh.XO4],[M0_L0_N,l0_bl_n,l0_bl_s,u0_bl_n_origin,u0_bl_s,M1_U0_S],[l0_w,l0_bl_w,u0_bl_w,u0_w],R_out,R_exit,bl_ratio,N_STATOR,Vertical_Line=Vertical_Line, alpha_exit=alpha_exit, arc=arc)
    else:
        print "ERROR: value of corner_flag is not supported."


    if corner_flag == 1:
        # M0
        M0_W = NGV_out.B7Xml#B9Xml#B9U3l#B9U3e#Line(NGV_out.B9,Mesh.X)
        M0_E = L0_W
        M0_N = NGV_out.XmlXm#Polyline([NGV_out.U3lU3,Line(Mesh.X,Mesh.Xm)],"",0.,1.)
        #B7B8 = Line(Mesh.XO4_BL,Mesh.XO3_BL)#Polyline([l0_bl_n],"",0.,l0_bl_n_ratio)
        M0_S = Line(Mesh.XO4_BL,Mesh.XO3_BL) #B7B8#Polyline([NGV_out.B9B1,B1B8],"",0.,1.)

        # M0_BL
        M0_BL_W = Line(Mesh.XO4,NGV_out.B7)
        M0_BL_E = L0_BL_W
        M0_BL_N = M0_S
        M0_BL_S = Line(Mesh.XO4,Mesh.XO3)

        # M1
        M1_W = NGV_out.YmlC6#D3lC6#D3eC6#D3C6#Arc(Mesh.Y,Mesh.YO2,Center)
        M1_E = U0_W
        M1_N = NGV_out.C6C4#T2B7
        M1_S = NGV_out.YmlYm#Polyline( [NGV_out.D3lD3,Line(Mesh.Y,Mesh.Ym)],"",0.,1.)

        # M1_BL  
        M1_BL_W = NGV_out.C6C5
        M1_BL_E = U0_BL_W
        M1_BL_N = NGV_out.C5C1
        M1_BL_S = M1_N    

        print "end of create line"
        #


    if gdata.dimensions == 2:
        if corner_flag == 1:
            #0 to 4
            print "creating L0_BL ..."
            L0_BL = Block2D(make_patch(L0_BL_N, L0_BL_E, L0_BL_S, L0_BL_W), nni=Nnoz0, nnj=Nbl,
                 cf_list=[No_30,BL_l,No_40,BL_l], fill_condition=initial, label="L0_BL")
            print "creating L0 ..."
            L0 = Block2D(make_patch(L0_N, L0_E, L0_S, L0_W), nni=Nnoz0, nnj=Nlow,
                 cf_list=[No_20,None,No_30,None], fill_condition=initial, label="L0")
            print "creating M0 ..."
            M0 = Block2D(make_patch(M0_N, M0_E, M0_S, M0_W), nni=Nmid, nnj=Nlow,
                 cf_list=[M0_M1,None,None,None], fill_condition=initial, label="M0")
            print "creating M0_BL ..."
            M0_BL = Block2D(make_patch(M0_BL_N, M0_BL_E, M0_BL_S, M0_BL_W), nni=Nmid, nnj=Nbl,
                   cf_list=[None,BL_l,None,BL_l], fill_condition=initial, label="M0_BL")   
            print "creating L1_BL ..."
            L1_BL = Block2D(make_patch(L1_BL_N, L1_BL_E, L1_BL_S, L1_BL_W), nni=Nnoz1, nnj=Nbl,
                  cf_list=[No_3,BL_l,No_4,BL_l],fill_condition=initial, label="L1_BL")

            #5 to 9
            print "creating L1 ..."
            L1 = Block2D(make_patch(L1_N, L1_E, L1_S, L1_W), nni=Nnoz1, nnj=Nlow,
                 cf_list=[No_2,L2_L1,No_3,None],fill_condition=initial, label="L1")
            print "creating L2 ..."
            L2 = Block2D(make_patch(L2_N, L2_E, L2_S, L2_W), nni=Nin, nnj=Nlow,
                 cf_list=[S_t2,None,L2_TOP,L2_L1],fill_condition=initial, label="L2")
            print "creating U0_BL ..."
            U0_BL = Block2D(make_patch(U0_BL_N, U0_BL_E, U0_BL_S, U0_BL_W), nni=Nnoz0, nnj=Nbl,
                 cf_list=[No_00,BL_r,No_10,BL_c],fill_condition=initial, label="U0_BL")
            print "creating U0 ..."
            U0 = Block2D(make_patch(U0_N, U0_E, U0_S, U0_W), nni=Nnoz0, nnj=Nupp,
                 cf_list=[No_10,None,No_20,M1_U0],fill_condition=initial, label="U0")
            print "creating M1 ..."
            M1 = Block2D(make_patch(M1_N, M1_E, M1_S, M1_W), nni=Nmid, nnj=Nupp,
                 cf_list=[None,M1_U0,M0_M1,S_out],fill_condition=initial, label="M1")  

            #10 to 15
            print "creating M1_BL ..."
            M1_BL = Block2D(make_patch(M1_BL_N, M1_BL_E, M1_BL_S, M1_BL_W), nni=Nmid, nnj=Nbl,
                 cf_list=[None,BL_c,None,BL_c],fill_condition=initial, label="M1_BL")  
            print "creating U1_BL ..."
            U1_BL = Block2D(make_patch(U1_BL_N, U1_BL_E, U1_BL_S, U1_BL_W), nni=Nnoz1, nnj=Nbl,
                 cf_list=[No_0,BL_r,No_1,BL_r],fill_condition=initial, label="U1_BL")
            print "creating U1 ..."
            U1 = Block2D(make_patch(U1_N, U1_E, U1_S, U1_W), nni=Nnoz1, nnj=Nupp,
                 cf_list=[No_1,None,No_2,None],fill_condition=initial, label="U1")
            print "creating U2 ..."
            U2 = Block2D(make_patch(U2_N, U2_E, U2_S, U2_W), nni=Nin, nnj=Nupp,
                 cf_list=[S_t0,None,S_t2,None],fill_condition=initial, label="U2")
            #print "creating TOP_BL ..."
            #TOP_BL = Block2D(make_patch(TOP_BL_N, TOP_BL_E, TOP_BL_S, TOP_BL_W), nni=Nbl, nnj=Ntip,
            #     cf_list=[BL_l,S_t1,BL_l,S_t1],fill_condition=initial, label="TOP_BL")
            print "creating TOPL_BL ..."
            TOPL_BL = Block2D(make_patch(TOPL_BL_N, TOPL_BL_E, TOPL_BL_S, TOPL_BL_W), nni=Nbl, nnj=Ntipl,
                 cf_list=[BL_l,S_t1,BL_l,S_t1],fill_condition=initial, label="TOPL_BL")
            print "creating TOPR_BL ..."
            TOPR_BL = Block2D(make_patch(TOPR_BL_N, TOPR_BL_E, TOPR_BL_S, TOPR_BL_W), nni=Nbl, nnj=Ntipr,
                 cf_list=[BL_l,S_t1,BL_l,S_t1],fill_condition=initial, label="TOPR_BL")

            #16 to 21
            print "creating TOPL ..."
            TOPL = Block2D(make_patch(TOPL_N, TOPL_E, TOPL_S, TOPL_W), nni=Nin, nnj=Ntipl,
                 cf_list=[L2_TOP,None,TOPR_TOPL,S_t1],fill_condition=initial, label="TOP0")
            print "creating TOPR ..."
            TOPR = Block2D(make_patch(TOPR_N, TOPR_E, TOPR_S, TOPR_W), nni=Nin, nnj=Ntipr,
                 cf_list=[TOPR_TOPL,None,S_t0,S_t1],fill_condition=initial, label="TOP0")
            print "creating E0 ..."
            E0 = Block2D(NGV_out.E0, nni=Nexit, nnj=Nmid,
                 cf_list=[None,None,None,None],fill_condition=initial, label="E0")
            print "creating E0_BL ..."
            E0_BL = Block2D(NGV_out.E0_BL, nni=Nbl, nnj=Nmid,
                 cf_list=[BL_r,None,BL_r,None],fill_condition=initial, label="E0_BL")
            print "creating E0_EX ..."
            E0_EX = Block2D(NGV_out.E0_EX, nni=Nout, nnj=Nmid,
                 cf_list=[None,None,None,None],fill_condition=initial, label="E0_EX")
            print "creating E1_BL ..."
            E1_BL = Block2D(NGV_out.E1_BL, nni=Nbl, nnj=Ntail0,
                 cf_list=[BL_r,None,BL_r,None],fill_condition=initial, label="E1_BL")

            #22 to 26
            print "creating E1 ..."
            E1 = Block2D(NGV_out.E1, nni=Nexit, nnj=Ntail0, 
                 cf_list=[None,None,None,None],fill_condition=initial, label="E1")
            print "creating E1_EX ..."
            E1_EX = Block2D(NGV_out.E1_EX, nni=Nout, nnj=Ntail0, 
                 cf_list=[None,None,None,None],fill_condition=initial, label="E1_EX")
            print "creating E2 ..."
            E2 = Block2D(NGV_out.E2, nni=Nexit, nnj=Ntail1, 
                 cf_list=[None,E2_m,E2_E3,E2_to],fill_condition=initial, label="E2")
            print "creating E2_BL ..."
            E2_BL = Block2D(NGV_out.E2_BL, nni=Nbl, nnj=Ntail1, 
                 cf_list=[BL_r,E2_tw,BL_r,E2_m],fill_condition=initial, label="E2_BL")
            print "creating E2_EX ..."
            E2_EX = Block2D(NGV_out.E2_EX, nni=Nout, nnj=Ntail1, 
                 cf_list=[None,E2_to,None,E2_o],fill_condition=initial, label="E2_EX")


            #27 to 29
            print "creating E3_BL ..."
            E3_BL = Block2D(NGV_out.E3_BL, nni=Nbl, nnj=Nmidd, 
                 cf_list=[BL_r,None,BL_c,None],fill_condition=initial, label="E3_BL")
            print "creating E3 ..."
            E3 = Block2D(NGV_out.E3, nni=Nexit, nnj=Nmidd, 
                 cf_list=[E2_E3,None,S_out,E3_Ex],fill_condition=initial, label="E3")
            print "creating E3_EX ..."
            E3_EX = Block2D(NGV_out.E3_EX, nni=Nout, nnj=Nmidd,
                 cf_list=[None,E3_Ex,None,E3_Ex],fill_condition=initial, label="E3_EX")
        elif corner_flag == 0:
        #0 to 4
            print "creating L0_BL ..."
            L0_BL = Block2D(make_patch(L0_BL_N, L0_BL_E, L0_BL_S, L0_BL_W), nni=Nnoz0, nnj=Nbl,
                        fill_condition=initial, label="L0_BL")
     
            L0 = Block2D(make_patch(L0_N, L0_E, L0_S, L0_W), nni=Nnoz0, nnj=Nlow,
                       fill_condition=initial, label="L0")

            L1_BL = Block2D(make_patch(L1_BL_N, L1_BL_E, L1_BL_S, L1_BL_W), nni=Nnoz1, nnj=Nbl,
                        fill_condition=initial, label="L1_BL")

            L1 = Block2D(make_patch(L1_N, L1_E, L1_S, L1_W), nni=Nnoz1, nnj=Nlow,
                        fill_condition=initial, label="L1")
     
            L2 = Block2D(make_patch(L2_N, L2_E, L2_S, L2_W), nni=Nin, nnj=Nlow,
                        fill_condition=initial, label="L2")
            ##
            U0_BL = Block2D(make_patch(U0_BL_N, U0_BL_E, U0_BL_S, U0_BL_W), nni=Nnoz0, nnj=Nbl,
                        fill_condition=initial, label="U0_BL")

     
            U0 = Block2D(make_patch(U0_N, U0_E, U0_S, U0_W), nni=Nnoz0, nnj=Nupp,
                       fill_condition=initial, label="U0")
     
            U1_BL = Block2D(make_patch(U1_BL_N, U1_BL_E, U1_BL_S, U1_BL_W), nni=Nnoz1, nnj=Nbl,
                        fill_condition=initial, label="U1_BL")

            U1 = Block2D(make_patch(U1_N, U1_E, U1_S, U1_W), nni=Nnoz1, nnj=Nupp,
                        fill_condition=initial, label="U1")

            U2 = Block2D(make_patch(U2_N, U2_E, U2_S, U2_W), nni=Nin, nnj=Nupp,
                        fill_condition=initial, label="U2")
            ##
            TOP_BL = Block2D(make_patch(TOP_BL_N, TOP_BL_E, TOP_BL_S, TOP_BL_W), nni=Nbl, nnj=Ntip,
                        fill_condition=initial, label="TOP_BL")

            TOP0 = Block2D(make_patch(TOP_N, TOP_E, TOP_S, TOP_W), nni=Nin, nnj=Ntip,
                        fill_condition=initial, label="TOP0")
            ##

            #E0 = Block2D(make_patch(data_E0[0],data_E0[1],data_E0[2],data_E0[3]), nni=Nexit, nnj=Nlow,
            #         fill_condition=initial, label="E0")
            E0 = Block2D(NGV_out.E0, nni=Nexit, nnj=Nlow,
                     fill_condition=initial, label="E0")

            E1 = Block2D(NGV_out.E1, nni=Nexit, nnj=Nbl,
                     fill_condition=initial, label="E1")

            E2 = Block2D(NGV_out.E2, nni=Nexit, nnj=Ntail,
                     fill_condition=initial, label="E2")

            E3 = Block2D(NGV_out.E3, nni=Nexit, nnj=Nbl,
                     fill_condition=initial, label="E3")

            E4 = Block2D(NGV_out.E4, nni=Nexit, nnj=Nupp,
                     fill_condition=initial, label="E4")
        else:
            print "Error, wrong value for corner_flag"

    if gdata.dimensions == 3:
        if corner_flag == 1:
            # block 0 to 3
            print "creating L0_BL ..."
            data = extruded_patch(make_patch(L0_BL_N, L0_BL_E, L0_BL_S, L0_BL_W), Vertical_Line)
            L0_BL = Block3D(data.Vol(), nni=Nnoz0, nnj=Nbl, nnk=Nz,
                 cf_list=[No_40,BL_l,No_30,BL_l,No_40,BL_l,No_30,BL_l,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="L0_BL")
            print "creating L0 ..."
            data = extruded_patch(make_patch(L0_N, L0_E, L0_S, L0_W), Vertical_Line)    
            L0 = Block3D( data.Vol(), nni=Nnoz0, nnj=Nlow, nnk=Nz,
                 cf_list=[No_30,None,No_20,None,No_30,None,No_20,None,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="L0") 
    

            ## blks 4 to 6
            print "creating L1_BL ..."
            data = extruded_patch(make_patch(L1_BL_N, L1_BL_E, L1_BL_S, L1_BL_W), Vertical_Line)
            L1_BL = Block3D(data.Vol(), nni=Nnoz1, nnj=Nbl, nnk=Nz,
                 cf_list=[None,BL_l,None,BL_l,None,BL_l,None,BL_l,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="L1_BL")
            print "creating L1 ..."
            data = extruded_patch(make_patch(L1_N, L1_E, L1_S, L1_W), Vertical_Line)
            L1 = Block3D(data.Vol(), nni=Nnoz1, nnj=Nlow, nnk=Nz,
                 cf_list=[None,L2_L1,No_1,None,None,L2_L1,No_1,None,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="L1")
            print "creating L2 ..."
            data = extruded_patch(make_patch(L2_N, L2_E, L2_S, L2_W), Vertical_Line)
            L2 = Block3D(data.Vol(), nni=Nin, nnj=Nlow, nnk=Nz,
                 cf_list=[L2_TOP,None,S_t2,L2_L1,L2_TOP,None,S_t2,L2_L1,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="L2")
            ## blks 7 to 8
            print "creating U0_BL ..."
            data = extruded_patch(make_patch(U0_BL_N, U0_BL_E, U0_BL_S, U0_BL_W), Vertical_Line)
            temp = data.Vol()
            #temp.rotate_about_zaxis(2.*np.pi/N_STATOR)
            U0_BL = Block3D(temp, nni=Nnoz0, nnj=Nbl, nnk=Nz,
                 cf_list=[No_00,BL_r,No_10,BL_c,No_00,BL_r,No_10,BL_c,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="U0_BL")
            print "creating U0 ..."
            data = extruded_patch(make_patch(U0_N, U0_E, U0_S, U0_W), Vertical_Line)
            temp = data.Vol()
            #temp.rotate_about_zaxis(2.*np.pi/N_STATOR)
            U0 = Block3D(temp, nni=Nnoz0, nnj=Nupp, nnk=Nz,
                 cf_list=[No_20,None,No_10,M1_U0,No_20,None,No_10,M1_U0,Z_d,Z_d,Z_d,Z_d],
                       fill_condition=initial, label="U0")

            ## blks 9 to 13
            print "creating U1_BL ..."
            data = extruded_patch(make_patch(U1_BL_N, U1_BL_E, U1_BL_S, U1_BL_W), Vertical_Line)
            U1_BL = Block3D(data.Vol(), nni=Nnoz1, nnj=Nbl, nnk=Nz,
                 cf_list=[No_1,BL_r,No_0,BL_r,No_1,BL_r,No_0,BL_r,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="U1_BL")
            print "creating U1 ..."
            data = extruded_patch(make_patch(U1_N, U1_E, U1_S, U1_W), Vertical_Line)
            U1 = Block3D(data.Vol(), nni=Nnoz1, nnj=Nupp, nnk=Nz,
                 cf_list=[No_0,None,No_1,None,No_0,None,No_1,None,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="U1")
            print "creating U2 ..."
            data = extruded_patch(make_patch(U2_N, U2_E, U2_S, U2_W), Vertical_Line)
            U2 = Block3D(data.Vol(), nni=Nin, nnj=Nupp, nnk=Nz,
                 cf_list=[S_t2,None,S_t0,None,S_t2,None,S_t0,None,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="U2")

            ## blks 14 to 15
            print "creating TOP_BL ..."
            data = extruded_patch(make_patch(TOP_BL_N, TOP_BL_E, TOP_BL_S, TOP_BL_W), Vertical_Line)
            TOP_BL = Block3D(data.Vol(), nni=Nbl, nnj=Ntip, nnk=Nz,
                 cf_list=[BL_l,S_t1,BL_l,S_t1,BL_l,S_t1,BL_l,S_t1,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="TOP_BL")
            print "creating TOP0 ..."
            data = extruded_patch(make_patch(TOP_N, TOP_E, TOP_S, TOP_W), Vertical_Line)
            TOP0 = Block3D(data.Vol(), nni=Nin, nnj=Ntip, nnk=Nz,
                 cf_list=[S_t0,None,L2_TOP,S_t1,S_t0,None,L2_TOP,S_t1,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="TOP0")
            print "creating M0 ..."
            data = extruded_patch(make_patch(M0_N, M0_E, M0_S, M0_W), Vertical_Line)
            M0 = Block3D(data.Vol(), nni=Nmid, nnj=Nlow, nnk=Nz,
                 cf_list=[None,None,M0_M1,None,None,None,M0_M1,None,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="M0")  
            print "creating M0_BL ..."
            data = extruded_patch(make_patch(M0_BL_N, M0_BL_E, M0_BL_S, M0_BL_W), Vertical_Line)
            M0_BL = Block3D(data.Vol(), nni=Nmid, nnj=Nbl, nnk=Nz,
                   cf_list=[None,BL_l,None,BL_l,None,BL_l,None,BL_l,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="M0_BL")    
            print "creating M1 ..."
            data = extruded_patch(make_patch(M1_N, M1_E, M1_S, M1_W), Vertical_Line)
            M1 = Block3D(data.Vol(), nni=Nmid, nnj=Nupp, nnk=Nz,
                 cf_list=[M0_M1,M1_U0,None,S_out,M0_M1,M1_U0,None,S_out,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="M1")  
            print "creating M1_BL ..."
            data = extruded_patch(make_patch(M1_BL_N, M1_BL_E, M1_BL_S, M1_BL_W), Vertical_Line)
            M1_BL = Block3D(data.Vol(), nni=Nmid, nnj=Nbl, nnk=Nz,
                 cf_list=[None,BL_c,None,BL_c,None,BL_c,None,BL_c,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="M1_BL")  
 
            ## blks 16 to 21
            print "creating E0 ..."
            data_E0 = NGV_out.eval_E0_3D()
            E0 = Block3D(data_E0, nni=Nexit, nnj=Nmid, nnk=Nz,
                 cf_list=[None,None,None,None,None,None,None,None,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E0")
            print "creating E0_BL ..."
            data_E0_BL = NGV_out.eval_E0_BL_3D()
            E0_BL = Block3D(data_E0_BL, nni=Nbl, nnj=Nmid, nnk=Nz,
                 cf_list=[BL_r,None,BL_r,None,BL_r,None,BL_r,None,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E0_BL")
            print "creating E0_EX ..."
            data_E0_EX = NGV_out.eval_E0_EX_3D()
            E0_EX = Block3D(data_E0_EX, nni=Nout, nnj=Nmid, nnk=Nz,
                 cf_list=[None,None,None,None,None,None,None,None,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E0_EX")
            print "creating E1_BL ..."
            data_E1_BL = NGV_out.eval_E1_BL_3D()
            E1_BL = Block3D(data_E1_BL, nni=Nbl, nnj=Ntail0, nnk=Nz,
                 cf_list=[BL_r,None,BL_r,None,BL_r,None,BL_r,None,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E1_BL")
            print "creating E1 ..."
            data_E1 = NGV_out.eval_E1_3D()
            E1 = Block3D(data_E1, nni=Nexit, nnj=Ntail0, nnk=Nz,
                 cf_list=[None,None,None,None,None,None,None,None,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E1")
            print "creating E1_EX ..."
            data_E1_EX = NGV_out.eval_E1_EX_3D()
            E1_EX = Block3D(data_E1_EX, nni=Nout, nnj=Ntail0, nnk=Nz,
                 cf_list=[None,None,None,None,None,None,None,None,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E1_EX")
            print "creating E2 ..."
            data_E2 = NGV_out.eval_E2_3D()
            E2 = Block3D(data_E2, nni=Nexit, nnj=Ntail1, nnk=Nz,
                 cf_list=[E2_E3,E2_m,None,E2_to,E2_E3,E2_m,None,E2_to,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E2")
            print "creating E2_BL ..."
            data_E2_BL = NGV_out.eval_E2_BL_3D()
            E2_BL = Block3D(data_E2_BL, nni=Nbl, nnj=Ntail1, nnk=Nz,
                 cf_list=[BL_r,E2_tw,BL_r,E2_m,BL_r,E2_tw,BL_r,E2_m,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E2_BL")
            print "creating E2_EX ..."
            data_E2_EX = NGV_out.eval_E2_EX_3D()
            E2_EX = Block3D(data_E2_EX, nni=Nout, nnj=Ntail1, nnk=Nz,
                 cf_list=[None,E2_to,None,E2_o,None,E2_to,None,E2_o,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E2_EX")
            print "creating E3_BL ..."
            data_E3_BL = NGV_out.eval_E3_BL_3D()
            E3_BL = Block3D(data_E3_BL, nni=Nbl, nnj=Nmidd, nnk=Nz,
                 cf_list=[BL_r,None,BL_c,None,BL_r,None,BL_c,None,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E3_BL")
            print "creating E3 ..."
            data_E3 = NGV_out.eval_E3_3D()
            E3 = Block3D(data_E3, nni=Nexit, nnj=Nmidd, nnk=Nz,
                 cf_list=[S_out,None,E2_E3,E3_Ex,S_out,None,E2_E3,E3_Ex,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E3")
            print "creating E3_EX ..."
            data_E3_EX = NGV_out.eval_E3_EX_3D()
            E3_EX = Block3D(data_E3_EX, nni=Nout, nnj=Nmidd, nnk=Nz,
                 cf_list=[None,E3_Ex,None,E3_Ex,None,E3_Ex,None,E3_Ex,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E3_EX")
            #print "creating E4 ..."       
            #data_E4 = NGV_out.eval_E4_3D()
            #E4 = Block3D(data_E4, nni=Nexit, nnj=Nmid, nnk=Nz,
            #     cf_list=[None,None,None,None,None,None,None,None,Z_d,Z_d,Z_d,Z_d],
            #         fill_condition=initial, label="E4")
            #print "creating E4_BL ..."       
            #data_E4_BL = NGV_out.eval_E4_BL_3D()
            #E4_BL = Block3D(data_E4_BL, nni=Nbl, nnj=Nmid, nnk=Nz,
            #     cf_list=[BL_r,None,BL_r,None,BL_r,None,BL_r,None,Z_d,Z_d,Z_d,Z_d],
            #         fill_condition=initial, label="E4_BL")
            #print "creating E5 ..."
            #data_E5 = NGV_out.eval_E5_3D()
            #E5 = Block3D(data_E5, nni=Nbl/2, nnj=Nmid, nnk=Nz,
            #     cf_list=[None,None,None,None,None,None,None,None,Z_d,Z_d,Z_d,Z_d],
            #         fill_condition=initial, label="E5")
        elif corner_flag ==  0:
            # block 0 to 3
            print "creating L0_BL ..."
            data = extruded_patch(make_patch(L0_BL_N, L0_BL_E, L0_BL_S, L0_BL_W), Vertical_Line)
            L0_BL = Block3D(data.Vol(), nni=Nnoz0, nnj=Nbl, nnk=Nz,
                 cf_list=[No_4,BL_l,No_3,BL_l,No_4,BL_l,No_3,BL_l,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="L0_BL")

            print "creating L0 ..."
            data = extruded_patch(make_patch(L0_N, L0_E, L0_S, L0_W), Vertical_Line)    
            L0 = Block3D( data.Vol(), nni=Nnoz0, nnj=Nlow, nnk=Nz,
                 cf_list=[No_3,None,No_2,None,No_3,None,No_2,None,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="L0") 
     
            print "creating L1_BL ..."
            data = extruded_patch(make_patch(L1_BL_N, L1_BL_E, L1_BL_S, L1_BL_W), Vertical_Line)
            L1_BL = Block3D(data.Vol(), nni=Nnoz1, nnj=Nbl, nnk=Nz,
                 cf_list=[None,BL_l,None,BL_l,None,BL_l,None,BL_l,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="L1_BL")

            print "creating L1 ..."
            data = extruded_patch(make_patch(L1_N, L1_E, L1_S, L1_W), Vertical_Line)
            L1 = Block3D(data.Vol(), nni=Nnoz1, nnj=Nlow, nnk=Nz,
                 cf_list=[None,None,None,None,None,None,None,None,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="L1")

            print "creating L2 ..."
            data = extruded_patch(make_patch(L2_N, L2_E, L2_S, L2_W), Vertical_Line)
            L2 = Block3D(data.Vol(), nni=Nin, nnj=Nlow, nnk=Nz,
                 cf_list=[S_t2,None,S_t0,None,S_t2,None,S_t0,None,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="L2")
            ## blks 5 to 9
            print "creating U0_BL ..."
            data = extruded_patch(make_patch(U0_BL_N, U0_BL_E, U0_BL_S, U0_BL_W), Vertical_Line)
            temp = data.Vol()
            temp.rotate_about_zaxis(2.*np.pi/N_STATOR)
            U0_BL = Block3D(temp, nni=Nnoz0, nnj=Nbl, nnk=Nz,
                 cf_list=[No_1,BL_r,No_0,BL_r,No_1,BL_r,No_0,BL_r,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="U0_BL")
     
            print "creating U0 ..."
            data = extruded_patch(make_patch(U0_N, U0_E, U0_S, U0_W), Vertical_Line)
            temp = data.Vol()
            temp.rotate_about_zaxis(2.*np.pi/N_STATOR)
            U0 = Block3D(temp, nni=Nnoz0, nnj=Nupp, nnk=Nz,
                 cf_list=[No_2,None,No_1,None,No_2,None,No_1,None,Z_d,Z_d,Z_d,Z_d],
                        fill_condition=initial, label="U0")
     
            print "creating U1_BL ..."
            data = extruded_patch(make_patch(U1_BL_N, U1_BL_E, U1_BL_S, U1_BL_W), Vertical_Line)
            U1_BL = Block3D(data.Vol(), nni=Nnoz1, nnj=Nbl, nnk=Nz,
                  cf_list=[No_1,BL_r,No_0,BL_r,No_1,BL_r,No_0,BL_r,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="U1_BL")
            print "creating U1 ..."
            data = extruded_patch(make_patch(U1_N, U1_E, U1_S, U1_W), Vertical_Line)
            U1 = Block3D(data.Vol(), nni=Nnoz1, nnj=Nupp, nnk=Nz,
                 cf_list=[None,None,No_1,None,None,None,No_1,None,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="U1")
            print "creating U2 ..."
            data = extruded_patch(make_patch(U2_N, U2_E, U2_S, U2_W), Vertical_Line)
            U2 = Block3D(data.Vol(), nni=Nin, nnj=Nupp, nnk=Nz,
                 cf_list=[S_t0,None,S_t2,None,S_t0,None,S_t2,None,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="U2")
            ## blks 14 to 15
            print "creating TOP_BL ..."
            data = extruded_patch(make_patch(TOP_BL_N, TOP_BL_E, TOP_BL_S, TOP_BL_W), Vertical_Line)
            TOP_BL = Block3D(data.Vol(), nni=Nbl, nnj=Ntip, nnk=Nz,
                  cf_list=[BL_l,S_t1,BL_l,S_t1,BL_l,S_t1,BL_l,S_t1,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="TOP_BL")
            print "creating TOP0 ..."
            data = extruded_patch(make_patch(TOP_N, TOP_E, TOP_S, TOP_W), Vertical_Line)
            TOP0 = Block3D(data.Vol(), nni=Nin, nnj=Ntip, nnk=Nz,
                 cf_list=[S_t2,None,S_t2,S_t1,S_t2,None,S_t2,S_t1,Z_d,Z_d,Z_d,Z_d],
                         fill_condition=initial, label="TOP0")

            print "creating E0 ..."
            data_E0 = NGV_out.eval_E0_3D()
            E0 = Block3D(data_E0, nni=Nexit, nnj=Nlow, nnk=Nz,
                  cf_list=[S_out,None,S_out,None,S_out,None,S_out,None,Z_d,Z_d,Z_d,Z_d],
                      fill_condition=initial, label="E0")
            print "creating E1 ..."
            data_E1 = NGV_out.eval_E1_3D()
            E1 = Block3D(data_E1, nni=Nexit, nnj=Nbl, nnk=Nz,
                 cf_list=[S_out,BL_l,S_out,BL_l,S_out,BL_l,S_out,BL_l,Z_d,Z_d,Z_d,Z_d],
                      fill_condition=initial, label="E1")
            print "creating E2 ..."
            data_E2 = NGV_out.eval_E2_3D()
            E2 = Block3D(data_E2, nni=Nexit, nnj=Ntail, nnk=Nz,
                 cf_list=[S_out,E2_tw,S_out,E2_to,S_out,E2_tw,S_out,E2_to,Z_d,Z_d,Z_d,Z_d],
                      fill_condition=initial, label="E2")
            print "creating E3 ..."
            data_E3 = NGV_out.eval_E3_3D()
            E3 = Block3D(data_E3, nni=Nexit, nnj=Nbl, nnk=Nz,
                 cf_list=[S_out,BL_r,S_out,BL_r,S_out,BL_r,S_out,BL_r,Z_d,Z_d,Z_d,Z_d],
                      fill_condition=initial, label="E3")
            print "creating E4 ..."
            data_E4 = NGV_out.eval_E4_3D()
            E4 = Block3D(data_E4, nni=Nexit, nnj=Nupp, nnk=Nz,
                 cf_list=[S_out,None,S_out,None,S_out,None,S_out,None,Z_d,Z_d,Z_d,Z_d],
                     fill_condition=initial, label="E4")
        else:
            print "Error, wrong value for corner_flag"

        if rear_cavity_flag == 1:

            def H_volume(data_G, R_exit, R_H, r,s,t):
                Coords = data_G.eval(0.,s,t) # get x,y,z coordinates along west edge of corresponding G-block
                #print Coords.x, Coords.y, Coords.z
                theta = np.arctan2(Coords.y,Coords.x)
                R = R_H + r * (R_exit - R_H)
                x = R * np.cos(theta)
                y = R * np.sin(theta)
                z = Coords.z
                return (x,y,z)
     
            Vertical_Line_F = Line(Node(0.,0.,Height0-Z_F), Node(0.,0.,Height0)) 
            NGV_F = NGV_outlet([Mesh.X,Mesh.XO1_BL,Mesh.XO1,Mesh.YO2,Mesh.YO2_BL,Mesh.Y],
                                [L0_N,L0_S,L0_BL_S,U0_BL_N,U0_N,U0_S],
                                [L0_W,L0_BL_W,U0_BL_W,U0_W],R_out,R_exit, 
                                Vertical_Line=Vertical_Line_F, alpha_exit=alpha_exit, arc=arc)
            Vertical_Line_G = Line(Node(0.,0.,Height0-Z_F-Z_G), Node(0.,0.,Height0-Z_F)) 
            NGV_G = NGV_outlet([Mesh.X,Mesh.XO1_BL,Mesh.XO1,Mesh.YO2,Mesh.YO2_BL,Mesh.Y],
                                [L0_N,L0_S,L0_BL_S,U0_BL_N,U0_N,U0_S],
                                [L0_W,L0_BL_W,U0_BL_W,U0_W],R_out,R_exit, 
                                Vertical_Line=Vertical_Line_G, alpha_exit=alpha_exit, arc=arc)
            ##
            cf2 = RobertsClusterFunction(1,1,1.1)
            cf3 = RobertsClusterFunction(0,1,1.1)
            cf4 = RobertsClusterFunction(0,1,1.1)
            data_F0 = NGV_F.eval_E0_3D()
            F0 = Block3D(data_F0, nni=Nexit, nnj=Nlow, nnk=nRC0,
                     cf_list=[None,None,None,None,None,None,None,None,cf2,cf2,cf2,cf2],
                     fill_condition=initial, label="F0")
            ##
            data_F1 = NGV_F.eval_E1_3D()
            F1 = Block3D(data_F1, nni=Nexit, nnj=Nbl, nnk=nRC0,
                     cf_list=[None,None,None,None,None,None,None,None,cf2,cf2,cf2,cf2],
                     fill_condition=initial, label="F1")
            ##
            data_F2 = NGV_F.eval_E2_3D()
            F2 = Block3D(data_F2, nni=Nexit, nnj=Ntail, nnk=nRC0,
                     cf_list=[None,None,None,None,None,None,None,None,cf2,cf2,cf2,cf2],
                     fill_condition=initial, label="F2")
            ##
            data_F3 = NGV_F.eval_E3_3D()
            F3 = Block3D(data_F3, nni=Nexit, nnj=Nbl, nnk=nRC0,
                     cf_list=[None,None,None,None,None,None,None,None,cf2,cf2,cf2,cf2],
                     fill_condition=initial, label="F3")
            ##
            data_F4 = NGV_F.eval_E4_3D()
            F4 = Block3D(data_F4, nni=Nexit, nnj=Nupp, nnk=nRC0,
                     cf_list=[None,None,None,None,None,None,None,None,cf2,cf2,cf2,cf2],
                     fill_condition=initial, label="F4")
            ##
            data_G0 = NGV_G.eval_E0_3D()
            G0 = Block3D(data_G0, nni=Nexit, nnj=Nlow, nnk=nRC1,
                     cf_list=[None,None,None,None,None,None,None,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="G0")
            ##
            data_G1 = NGV_G.eval_E1_3D()
            G1 = Block3D(data_G1, nni=Nexit, nnj=Nbl, nnk=nRC1,
                     cf_list=[None,None,None,None,None,None,None,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="G1")
            ##
            data_G2 = NGV_G.eval_E2_3D()
            G2 = Block3D(data_G2, nni=Nexit, nnj=Ntail, nnk=nRC1,
                     cf_list=[None,None,None,None,None,None,None,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="G2")
            ##
            data_G3 = NGV_G.eval_E3_3D()
            G3 = Block3D(data_G3, nni=Nexit, nnj=Nbl, nnk=nRC1,
                     cf_list=[None,None,None,None,None,None,None,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="G3")
            ##
            data_G4 = NGV_G.eval_E4_3D()
            G4 = Block3D(data_G4, nni=Nexit, nnj=Nupp, nnk=nRC1,
                     cf_list=[None,None,None,None,None,None,None,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="G4")
            ##
            pyfunction_H0 = lambda r,s,t: H_volume(data_G0, R_exit, R_H, r,s,t)
            H0 = Block3D(PyFunctionVolume(pyfunction_H0), nni=nRC2, nnj=Nlow, nnk=nRC1,
                     cf_list=[cf4,None,cf4,None,cf4,None,cf4,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="H0")
            ##
            pyfunction_H1 = lambda r,s,t: H_volume(data_G1, R_exit, R_H, r,s,t)
            H1 = Block3D(PyFunctionVolume(pyfunction_H1), nni=nRC2, nnj=Nbl, nnk=nRC1,
                     cf_list=[cf4,None,cf4,None,cf4,None,cf4,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="H1")
            ##
            pyfunction_H2 = lambda r,s,t: H_volume(data_G2, R_exit, R_H, r,s,t)
            H2 = Block3D(PyFunctionVolume(pyfunction_H2), nni=nRC2, nnj=Ntail, nnk=nRC1,
                     cf_list=[cf4,None,cf4,None,cf4,None,cf4,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="H2")
            ##
            pyfunction_H3 = lambda r,s,t: H_volume(data_G3, R_exit, R_H, r,s,t)
            H3 = Block3D(PyFunctionVolume(pyfunction_H3), nni=nRC2, nnj=Nbl, nnk=nRC1,
                     cf_list=[cf4,None,cf4,None,cf4,None,cf4,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="H3")
            ##
            pyfunction_H4 = lambda r,s,t: H_volume(data_G4, R_exit, R_H, r,s,t)
            H4 = Block3D(PyFunctionVolume(pyfunction_H4), nni=nRC2, nnj=Nupp, nnk=nRC1,
                     cf_list=[cf4,None,cf4,None,cf4,None,cf4,None,cf3,cf3,cf3,cf3],
                     fill_condition=initial, label="H4")
        
            """
            Z_F = 2.e-3
            Z_G = 2.e-3
            R_H = R_exit - 4.e-3
            nRC0 = 10
            nRC1 = 10
            nRC2 = 20
            """
        
    # link blocks
    identify_block_connections()

    if solver_type is "Eilmer":
        # Connect Periodic Blocks
        theta = 2.*np.pi/N_STATOR
        connect_blocks_3D(E4,E0,[(1,2),(5,6),(4,7),(0,3)],
                reorient_vector_quantities=True, 
                nA=[0.,1.,0.],t1A=[1.,0.,0.],
                nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
                check_corner_locations=False)

        connect_blocks_3D(U0,L0,[(1,2),(5,6),(4,7),(0,3)],
                reorient_vector_quantities=True, 
                nA=[0.,1.,0.],t1A=[1.,0.,0.],
                nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
                check_corner_locations=False)

        connect_blocks_3D(U1,L1,[(1,2),(5,6),(4,7),(0,3)],
                reorient_vector_quantities=True, 
                nA=[0.,1.,0.],t1A=[1.,0.,0.],
                nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
                check_corner_locations=False)

        connect_blocks_3D(U2,L2,[(1,2),(5,6),(4,7),(0,3)],
                reorient_vector_quantities=True, 
                nA=[0.,1.,0.],t1A=[1.,0.,0.],
                nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
                check_corner_locations=False)

        if rear_cavity_flag == 1:
            connect_blocks_3D(F4,F0,[(1,2),(5,6),(4,7),(0,3)],
                    reorient_vector_quantities=True, 
                    nA=[0.,1.,0.],t1A=[1.,0.,0.],
                    nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
                    check_corner_locations=False)
            connect_blocks_3D(G4,G0,[(1,2),(5,6),(4,7),(0,3)],
                    reorient_vector_quantities=True, 
                    nA=[0.,1.,0.],t1A=[1.,0.,0.],
                    nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
                    check_corner_locations=False)
            connect_blocks_3D(H4,H0,[(1,2),(5,6),(4,7),(0,3)],
                    reorient_vector_quantities=True, 
                    nA=[0.,1.,0.],t1A=[1.,0.,0.],
                    nB=[-np.sin(theta),np.cos(theta),0.],t1B=[np.cos(theta),np.sin(theta),0.],
                    check_corner_locations=False)

    if solver_type is "OF":
        if corner_flag == 1:
            # set correct labels for Periodic Boundaries.
            E0_EX.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            L0.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            L1.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            L2.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            M0.bc_list[NORTH] =  SlipWallBC(label='OF_symmetry_00')
            
            E3_EX.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_01')
            U0.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_01')
            U1.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_01')
            U2.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_01')
            M1.bc_list[SOUTH] =  SlipWallBC(label='OF_symmetry_01') 
        elif corner_flag == 0:
            # set correct labels for Periodic Boundaries.
            E0.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            L0.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            L1.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            L2.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00') 
             
            E4.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_01')
            U0.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_01')
            U1.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_01')
            U2.bc_list[SOUTH] = SlipWallBC(label='OF_symmetry_01')    

        else:
            print "Error"   

        if rear_cavity_flag == 1:
            F0.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            G0.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            H0.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')            

            F4.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            G4.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')
            H4.bc_list[NORTH] = SlipWallBC(label='OF_symmetry_00')   
            
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

    L2.bc_list[EAST] = SubsonicInBC(stagnation,label='OF_inlet_00',
                        direction_type='radial',direction_alpha=inlet_flow_direction)
    TOPL.bc_list[EAST] = SubsonicInBC(stagnation,label='OF_inlet_00',
                        direction_type='radial',direction_alpha=inlet_flow_direction)
    TOPR.bc_list[EAST] = SubsonicInBC(stagnation,label='OF_inlet_00',
                        direction_type='radial',direction_alpha=inlet_flow_direction)
    U2.bc_list[EAST] = SubsonicInBC(stagnation,label='OF_inlet_00',
                        direction_type='radial',direction_alpha=inlet_flow_direction)

    if corner_flag == 1:
        E0_EX.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00')
        E1_EX.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00')
        E2_EX.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00') 
        E3_EX.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00')


        L0_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
        L1_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
        M0_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
        E0_BL.bc_list[EAST] = AdiabaticBC(label='OF_wall_00')
        U0_BL.bc_list[NORTH] = AdiabaticBC(label='OF_wall_00')
        U1_BL.bc_list[NORTH] = AdiabaticBC(label='OF_wall_00')
        M1_BL.bc_list[NORTH] = AdiabaticBC(label='OF_wall_00')
        TOPL_BL.bc_list[WEST] = AdiabaticBC(label='OF_wall_00')
        TOPR_BL.bc_list[WEST] = AdiabaticBC(label='OF_wall_00')

        E1_BL.bc_list[EAST] = AdiabaticBC(label='OF_wall_01')
        E2_BL.bc_list[EAST] = AdiabaticBC(label='OF_wall_01')
        E3_BL.bc_list[EAST] = AdiabaticBC(label='OF_wall_01')
    elif corner_flag == 0:
        #E0.bc_list[WEST] = ExtrapolateOutBC(label='OF_outlet_00')
        #E1.bc_list[WEST] = ExtrapolateOutBC(label='OF_outlet_00')
        #E2.bc_list[WEST] = ExtrapolateOutBC(label='OF_outlet_00')
        #E3.bc_list[WEST] = ExtrapolateOutBC(label='OF_outlet_00')
        #E4.bc_list[WEST] = ExtrapolateOutBC(label='OF_outlet_00')

        E0.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00')
        E1.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00')
        E2.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00')
        E3.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00')
        E4.bc_list[WEST] = FixedPOutBC(Pout, label='OF_outlet_00')

        #L0_BL.bc_list[SOUTH] = SlipWallBC(label='OF_wall_00')
        #L1_BL.bc_list[SOUTH] = SlipWallBC(label='OF_wall_00')
        #U0_BL.bc_list[NORTH] = SlipWallBC(label='OF_wall_00')
        #U1_BL.bc_list[NORTH] = SlipWallBC(label='OF_wall_00')
        #TOP_BL.bc_list[WEST] = SlipWallBC(label='OF_wall_00')
        #E2.bc_list[EAST] = SlipWallBC(label='OF_wall_01')

        L0_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
        L1_BL.bc_list[SOUTH] = AdiabaticBC(label='OF_wall_00')
        U0_BL.bc_list[NORTH] = AdiabaticBC(label='OF_wall_00')
        U1_BL.bc_list[NORTH] = AdiabaticBC(label='OF_wall_00')
        TOP_BL.bc_list[WEST] = AdiabaticBC(label='OF_wall_00')
        E2.bc_list[EAST] = AdiabaticBC(label='OF_wall_01')

    else:
        print "Error"

    if gdata.dimensions == 3:
        if corner_flag == 1:
            E0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E0_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E0_EX.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E1.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E1_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E1_EX.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E2.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E2_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E2_EX.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E3.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E3_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E3_EX.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')

            if rear_cavity_flag == 0:
                E0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E0_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E0_EX.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E1_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E1_EX.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E2.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E2_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E2_EX.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E3.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E3_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E3_EX.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            L0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            L1.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            L2.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')   
            L0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            L1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            L2.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            L0_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            L1_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            L0_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            L1_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            U0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            U1.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            U2.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')   
            U0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            U1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            U2.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            U0_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            U1_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            U0_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            U1_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            TOP0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            TOP0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            TOP_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            TOP_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            
            M0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            M0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            M0_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            M0_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            M1.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            M1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            M1_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            M1_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')      
        elif corner_flag == 0:
            E0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E1.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E2.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E3.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            E4.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            if rear_cavity_flag == 0:
                E0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E2.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E3.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
                E4.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            L0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            L1.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            L2.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')   
            L0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            L1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            L2.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            L0_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            L1_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            L0_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            L1_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            U0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            U1.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            U2.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')   
            U0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            U1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            U2.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            U0_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            U1_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            U0_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
            U1_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            TOP0.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            TOP0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')

            TOP_BL.bc_list[TOP] = AdiabaticBC(label='OF_wall_02')
            TOP_BL.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_03')
        else:
            print "Error"

 
        if rear_cavity_flag == 1:
            F0.bc_list[WEST] = AdiabaticBC(label='OF_wall_04')
            F1.bc_list[WEST] = AdiabaticBC(label='OF_wall_04')
            F2.bc_list[WEST] = AdiabaticBC(label='OF_wall_04')
            F3.bc_list[WEST] = AdiabaticBC(label='OF_wall_04')
            F4.bc_list[WEST] = AdiabaticBC(label='OF_wall_04')
            ##
            H0.bc_list[TOP] = AdiabaticBC(label='OF_wall_04')
            H1.bc_list[TOP] = AdiabaticBC(label='OF_wall_04')
            H2.bc_list[TOP] = AdiabaticBC(label='OF_wall_04')
            H3.bc_list[TOP] = AdiabaticBC(label='OF_wall_04')
            H4.bc_list[TOP] = AdiabaticBC(label='OF_wall_04')
            ##
            F0.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            F1.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            F2.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            F3.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            F4.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            ##
            G0.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            G1.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            G2.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            G3.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            G4.bc_list[EAST] = AdiabaticBC(label='OF_wall_05')
            ##
            G0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            G1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            G2.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            G3.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            G4.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            ##
            H0.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            H1.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            H2.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            H3.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            H4.bc_list[BOTTOM] = AdiabaticBC(label='OF_wall_05')
            ##
            H0.bc_list[WEST] = FixedPOutBC(Pout_RC, label='OF_outlet_01')
            H1.bc_list[WEST] = FixedPOutBC(Pout_RC, label='OF_outlet_01')
            H2.bc_list[WEST] = FixedPOutBC(Pout_RC, label='OF_outlet_01')
            H3.bc_list[WEST] = FixedPOutBC(Pout_RC, label='OF_outlet_01')
            H4.bc_list[WEST] = FixedPOutBC(Pout_RC, label='OF_outlet_01')

    Coords = Mesh.Xi
    r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
    T1_in = t
    Coords = Mesh.Yi
    r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
    T0_in = t

    if corner_flag == 1:
        Coords = NGV_out.E0_EX.eval(0.,1.)
        r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
        T1_out = t
        Coords = NGV_out.E3_EX.eval(0.,0.)
        r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
        T0_out = t
    elif corner_flag == 0:
        Coords = NGV_out.E0.eval(0.,1.)
        r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
        T1_out = t
        Coords = NGV_out.E4.eval(0.,0.)
        r,t,z = cart2polar(Coords.x,Coords.y,Coords.z)
        T0_out = t

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

    #sketch.prefer_bc_labels_on_faces() # required to allow grouping of boundaries by e3prepToFoam.py

    if gdata.dimensions == 2:
       print "-------------------------"
        # This is to make a nice *.svg file of the 2-D projection of the mesh
       # sketch.xaxis(-0.05, 0.05, 0.02, 0.0)
       # sketch.yaxis(0.05, 0.05, 0.02, 0.0)
       # sketch.window(0., -0.025, 0.05, 0.025, 0.01, 0.01, 0.2, 0.2) 
       sketch.window(xmin=0.065, ymin=0.002, xmax=0.064, ymax=0.0025)

