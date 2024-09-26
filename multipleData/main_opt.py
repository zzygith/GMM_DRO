###### x1+x2+x1*xi
########################## x variable
######################################### for loop
######################################################## pi range
from pyomo.environ import *
from scipy.optimize import linprog
import copy
import numpy as np

# 定义系数矩阵
A_eq = [
    [1, 1, 1, 0, 0, 0, 0, 0, 0],  # pi11 + pi12 + pi13 = 0.6
    [0, 0, 0, 1, 1, 1, 0, 0, 0],  # pi21 + pi22 + pi23 = 0.25
    [0, 0, 0, 0, 0, 0, 1, 1, 1],  # pi31 + pi32 + pi33 = 0.15
    [1, 0, 0, 1, 0, 0, 1, 0, 0],  # pi11 + pi21 + pi31 = 0.6
    [0, 1, 0, 0, 1, 0, 0, 1, 0],  # pi12 + pi22 + pi32 = 0.25
    [0, 0, 1, 0, 0, 1, 0, 0, 1]   # pi13 + pi23 + pi33 = 0.15
]

#b_eq = [0.6, 0.25, 0.15, 0.6, 0.25, 0.15]
b_eq = [0.126, 0.622, 0.252, 0.126, 0.622, 0.252]

# 目标函数的系数列表（用于最大化和最小化每个变量）
c_list = [
    [1, 0, 0, 0, 0, 0, 0, 0, 0],  # Maximize or minimize a11
    [0, 1, 0, 0, 0, 0, 0, 0, 0],  # Maximize or minimize a12
    [0, 0, 1, 0, 0, 0, 0, 0, 0],  # Maximize or minimize a13
    [0, 0, 0, 1, 0, 0, 0, 0, 0],  # Maximize or minimize a21
    [0, 0, 0, 0, 1, 0, 0, 0, 0],  # Maximize or minimize a22
    [0, 0, 0, 0, 0, 1, 0, 0, 0],  # Maximize or minimize a23
    [0, 0, 0, 0, 0, 0, 1, 0, 0],  # Maximize or minimize a31
    [0, 0, 0, 0, 0, 0, 0, 1, 0],  # Maximize or minimize a32
    [0, 0, 0, 0, 0, 0, 0, 0, 1]   # Maximize or minimize a33
]

# 结果字典
resultsRange = {}

# 对每个变量进行最大化和最小化
for i, c in enumerate(c_list):
    # 最小化
    res_min = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0, None))
    # 最大化
    res_max = linprog([-x for x in c], A_eq=A_eq, b_eq=b_eq, bounds=(0, None))

    # 存储结果
    row = i // 3 + 1
    col = i % 3 + 1
    resultsRange[f'pi{row}{col}'] = (res_min.fun, -res_max.fun)


boundList={}
#distanceThreshold=100.0
distanceThreshold=10.0
#distanceThreshold=1.0
#distanceThreshold=0.01

# bUP1=[1.67, 1320.41, 5077.55,   1269.39, 4, 1363.05,  4895.23,1320.41,6.66668]
# bLOW1=[0.0,  1045.71, 4368.28,   1092.07,0.0,1008.41,  4540.60,1045.71,0.0]
# boundList[1.0]=[bUP1,bLOW1]

bUP1=[7.93652, 4956.71, 1338.28,    5177.25, 1.60773, 1338.37,   1399.33, 1285.92, 3.96826]
bLOW1=[0.0,    4606.06, 1062.66,    4398.18, 0.0,     1062.74,   1009.54, 1110.47, 0.0]
boundList[1.0]=[bUP1,bLOW1]

bUP1=[79.3651, 5350.27, 1671.98,    6090.97, 16.0772, 1672.08,   1892.18, 1490.06, 39.6825]
bLOW1=[0.0,    4241.43, 800.385,    3627.32, 0.0,     800.454,   659.553, 935.263, 0.0]
boundList[10.0]=[bUP1,bLOW1]

boundUpperList=boundList[distanceThreshold][0]
boundLowerList=boundList[distanceThreshold][1]

boundUpperList_copy = copy.deepcopy(boundUpperList)
boundLowerList_copy = copy.deepcopy(boundLowerList)


def primalFunc(xVectorSection,resultsRange,boundUpperListF,boundLowerListF):
    ##################################################################### d11 upper bound, pi parameter, UB
    model = ConcreteModel()
    pi11L=resultsRange['pi11'][0]
    pi12L=resultsRange['pi12'][0]
    pi13L=resultsRange['pi13'][0]
    pi21L=resultsRange['pi21'][0]
    pi22L=resultsRange['pi22'][0]
    pi23L=resultsRange['pi23'][0]
    pi31L=resultsRange['pi31'][0]
    pi32L=resultsRange['pi32'][0]
    pi33L=resultsRange['pi33'][0]

    pi11U=resultsRange['pi11'][1]
    pi12U=resultsRange['pi12'][1]
    pi13U=resultsRange['pi13'][1]
    pi21U=resultsRange['pi21'][1]
    pi22U=resultsRange['pi22'][1]
    pi23U=resultsRange['pi23'][1]
    pi31U=resultsRange['pi31'][1]
    pi32U=resultsRange['pi32'][1]
    pi33U=resultsRange['pi33'][1]

    d11U=boundUpperListF[0]
    d12U=boundUpperListF[1]
    d13U=boundUpperListF[2]
    d21U=boundUpperListF[3]
    d22U=boundUpperListF[4]
    d23U=boundUpperListF[5]
    d31U=boundUpperListF[6]
    d32U=boundUpperListF[7]
    d33U=boundUpperListF[8]

    d11L=boundLowerListF[0]
    d12L=boundLowerListF[1]
    d13L=boundLowerListF[2]
    d21L=boundLowerListF[3]
    d22L=boundLowerListF[4]
    d23L=boundLowerListF[5]
    d31L=boundLowerListF[6]
    d32L=boundLowerListF[7]
    d33L=boundLowerListF[8]

    model.pi11=Var(within=NonNegativeReals, initialize=wLH1)
    model.pi12=Var(within=NonNegativeReals, initialize=0.0)
    model.pi13=Var(within=NonNegativeReals, initialize=0.0)
    model.pi21=Var(within=NonNegativeReals, initialize=0.0)
    model.pi22=Var(within=NonNegativeReals, initialize=wLH2)
    model.pi23=Var(within=NonNegativeReals, initialize=0.0)
    model.pi31=Var(within=NonNegativeReals, initialize=0.0)
    model.pi32=Var(within=NonNegativeReals, initialize=0.0)
    model.pi33=Var(within=NonNegativeReals, initialize=wLH3)

    pi11=model.pi11
    pi12=model.pi12
    pi13=model.pi13
    pi21=model.pi21
    pi22=model.pi22
    pi23=model.pi23
    pi31=model.pi31
    pi32=model.pi32
    pi33=model.pi33


    model.ec1 = Constraint(expr = pi11+pi12+pi13==wUH1)
    model.ec2 = Constraint(expr = pi21+pi22+pi23==wUH2)
    model.ec3 = Constraint(expr = pi31+pi32+pi33==wUH3)
    model.ec4 = Constraint(expr = pi11+pi21+pi31==wLH1)
    model.ec5 = Constraint(expr = pi12+pi22+pi32==wLH2)
    model.ec6 = Constraint(expr = pi13+pi23+pi33==wLH3)


    # model.miuLH1=Var( within=Reals, initialize=miuUH1)
    # model.alphaLH1=Var(within=NonNegativeReals, initialize=alphaUH1)
    # model.miuLH2=Var( within=Reals, initialize=miuUH2)
    # model.alphaLH2=Var(within=NonNegativeReals, initialize=alphaUH2)
    # model.miuLH3=Var( within=Reals, initialize=miuUH3)
    # model.alphaLH3=Var(within=NonNegativeReals, initialize=alphaUH3)

    # miuLH1=model.miuLH1
    # alphaLH1=model.alphaLH1
    # miuLH2=model.miuLH2
    # alphaLH2=model.alphaLH2
    # miuLH3=model.miuLH3
    # alphaLH3=model.alphaLH3

    def miuLH1_init_rule(model, i):
        return miuUH1[i]
    def miuLH2_init_rule(model, i):
        return miuUH2[i]
    def miuLH3_init_rule(model, i):
        return miuUH3[i]
    def sigmaLH1_init_rule(model, i):
        return sigmaUH1[i]
    def sigmaLH2_init_rule(model, i):
        return sigmaUH2[i]
    def sigmaLH3_init_rule(model, i):
        return sigmaUH3[i]

    # 创建变量数组
    model.indices = RangeSet(0, 9)
    model.miuLH1 = Var(model.indices, within=Reals, initialize=miuLH1_init_rule)
    model.miuLH2 = Var(model.indices, within=Reals, initialize=miuLH2_init_rule)
    model.miuLH3 = Var(model.indices, within=Reals, initialize=miuLH3_init_rule)
    model.sigmaLH1=Var(model.indices, within=NonNegativeReals, initialize=sigmaLH1_init_rule)
    model.sigmaLH2=Var(model.indices, within=NonNegativeReals, initialize=sigmaLH2_init_rule)
    model.sigmaLH3=Var(model.indices, within=NonNegativeReals, initialize=sigmaLH3_init_rule)

    miuLH1=model.miuLH1
    sigmaLH1=model.sigmaLH1
    miuLH2=model.miuLH2
    sigmaLH2=model.sigmaLH2
    miuLH3=model.miuLH3
    sigmaLH3=model.sigmaLH3

    model.d11=Var(within=NonNegativeReals, initialize=0.0)
    model.d12=Var(within=NonNegativeReals, initialize=2.0)
    model.d13=Var(within=NonNegativeReals, initialize=0.0)
    model.d21=Var(within=NonNegativeReals, initialize=2.0)
    model.d22=Var(within=NonNegativeReals, initialize=0.0)
    model.d23=Var(within=NonNegativeReals, initialize=0.0)
    model.d31=Var(within=NonNegativeReals, initialize=0.0)
    model.d32=Var(within=NonNegativeReals, initialize=0.0)
    model.d33=Var(within=NonNegativeReals, initialize=0.0)

    d11=model.d11
    d12=model.d12
    d13=model.d13
    d21=model.d21
    d22=model.d22
    d23=model.d23
    d31=model.d31
    d32=model.d32
    d33=model.d33

    model.w11=Var(within=NonNegativeReals, initialize=0.0)
    model.w12=Var(within=NonNegativeReals, initialize=0.0)
    model.w13=Var(within=NonNegativeReals, initialize=0.0)
    model.w21=Var(within=NonNegativeReals, initialize=0.0)
    model.w22=Var(within=NonNegativeReals, initialize=0.0)
    model.w23=Var(within=NonNegativeReals, initialize=0.0)
    model.w31=Var(within=NonNegativeReals, initialize=0.0)
    model.w32=Var(within=NonNegativeReals, initialize=0.0)
    model.w33=Var(within=NonNegativeReals, initialize=0.0)

    w11=model.w11
    w12=model.w12
    w13=model.w13
    w21=model.w21
    w22=model.w22
    w23=model.w23
    w31=model.w31
    w32=model.w32
    w33=model.w33


    # model.ine1_11=Constraint(expr=d11-((miuUH1-miuLH1)**2+alphaUH1**2+alphaLH1**2-2*alphaUH1*alphaLH1)>=0)
    # model.ine1_12=Constraint(expr=d12-((miuUH1-miuLH2)**2+alphaUH1**2+alphaLH2**2-2*alphaUH1*alphaLH2)>=0)
    # model.ine1_13=Constraint(expr=d13-((miuUH1-miuLH3)**2+alphaUH1**2+alphaLH3**2-2*alphaUH1*alphaLH3)>=0)
    # model.ine1_21=Constraint(expr=d21-((miuUH2-miuLH1)**2+alphaUH2**2+alphaLH1**2-2*alphaUH2*alphaLH1)>=0)
    # model.ine1_22=Constraint(expr=d22-((miuUH2-miuLH2)**2+alphaUH2**2+alphaLH2**2-2*alphaUH2*alphaLH2)>=0)
    # model.ine1_23=Constraint(expr=d23-((miuUH2-miuLH3)**2+alphaUH2**2+alphaLH3**2-2*alphaUH2*alphaLH3)>=0)
    # model.ine1_31=Constraint(expr=d31-((miuUH3-miuLH1)**2+alphaUH3**2+alphaLH1**2-2*alphaUH3*alphaLH1)>=0)
    # model.ine1_32=Constraint(expr=d32-((miuUH3-miuLH2)**2+alphaUH3**2+alphaLH2**2-2*alphaUH3*alphaLH2)>=0)
    # model.ine1_33=Constraint(expr=d33-((miuUH3-miuLH3)**2+alphaUH3**2+alphaLH3**2-2*alphaUH3*alphaLH3)>=0)
    model.ine1_11=Constraint(expr=d11-( sum([ (miuUH1[i]-miuLH1[i])**2 for i in model.indices])+ sum( [sigmaUH1[i]+sigmaLH1[i]-2*(sigmaUH1[i]*sigmaLH1[i])**0.5 for i in model.indices] ) )>=0)
    model.ine1_12=Constraint(expr=d12-( sum([ (miuUH1[i]-miuLH2[i])**2 for i in model.indices])+ sum( [sigmaUH1[i]+sigmaLH2[i]-2*(sigmaUH1[i]*sigmaLH2[i])**0.5 for i in model.indices] ) )>=0)
    model.ine1_13=Constraint(expr=d13-( sum([ (miuUH1[i]-miuLH3[i])**2 for i in model.indices])+ sum( [sigmaUH1[i]+sigmaLH3[i]-2*(sigmaUH1[i]*sigmaLH3[i])**0.5 for i in model.indices] ) )>=0)
    model.ine1_21=Constraint(expr=d21-( sum([ (miuUH2[i]-miuLH1[i])**2 for i in model.indices])+ sum( [sigmaUH2[i]+sigmaLH1[i]-2*(sigmaUH2[i]*sigmaLH1[i])**0.5 for i in model.indices] ) )>=0)
    model.ine1_22=Constraint(expr=d22-( sum([ (miuUH2[i]-miuLH2[i])**2 for i in model.indices])+ sum( [sigmaUH2[i]+sigmaLH2[i]-2*(sigmaUH2[i]*sigmaLH2[i])**0.5 for i in model.indices] ) )>=0)
    model.ine1_23=Constraint(expr=d23-( sum([ (miuUH2[i]-miuLH3[i])**2 for i in model.indices])+ sum( [sigmaUH2[i]+sigmaLH3[i]-2*(sigmaUH2[i]*sigmaLH3[i])**0.5 for i in model.indices] ) )>=0)
    model.ine1_31=Constraint(expr=d31-( sum([ (miuUH3[i]-miuLH1[i])**2 for i in model.indices])+ sum( [sigmaUH3[i]+sigmaLH1[i]-2*(sigmaUH3[i]*sigmaLH1[i])**0.5 for i in model.indices] ) )>=0)
    model.ine1_32=Constraint(expr=d32-( sum([ (miuUH3[i]-miuLH2[i])**2 for i in model.indices])+ sum( [sigmaUH3[i]+sigmaLH2[i]-2*(sigmaUH3[i]*sigmaLH2[i])**0.5 for i in model.indices] ) )>=0)
    model.ine1_33=Constraint(expr=d33-( sum([ (miuUH3[i]-miuLH3[i])**2 for i in model.indices])+ sum( [sigmaUH3[i]+sigmaLH3[i]-2*(sigmaUH3[i]*sigmaLH3[i])**0.5 for i in model.indices] ) )>=0)

    model.ineWLG11=Constraint(expr=w11>=pi11L*d11+pi11*d11L-pi11L*d11L)
    model.ineWLG12=Constraint(expr=w12>=pi12L*d12+pi12*d12L-pi12L*d12L)
    model.ineWLG13=Constraint(expr=w13>=pi13L*d13+pi13*d13L-pi13L*d13L)
    model.ineWLG21=Constraint(expr=w21>=pi21L*d21+pi21*d21L-pi21L*d21L)
    model.ineWLG22=Constraint(expr=w22>=pi22L*d22+pi22*d22L-pi22L*d22L)
    model.ineWLG23=Constraint(expr=w23>=pi23L*d23+pi23*d23L-pi23L*d23L)
    model.ineWLG31=Constraint(expr=w31>=pi31L*d31+pi31*d31L-pi31L*d31L)
    model.ineWLG32=Constraint(expr=w32>=pi32L*d32+pi32*d32L-pi32L*d32L)
    model.ineWLG33=Constraint(expr=w33>=pi33L*d33+pi33*d33L-pi33L*d33L)

    model.ineWUG11=Constraint(expr=w11>=pi11U*d11+pi11*d11U-pi11U*d11U)
    model.ineWUG12=Constraint(expr=w12>=pi12U*d12+pi12*d12U-pi12U*d12U)
    model.ineWUG13=Constraint(expr=w13>=pi13U*d13+pi13*d13U-pi13U*d13U)
    model.ineWUG21=Constraint(expr=w21>=pi21U*d21+pi21*d21U-pi21U*d21U)
    model.ineWUG22=Constraint(expr=w22>=pi22U*d22+pi22*d22U-pi22U*d22U)
    model.ineWUG23=Constraint(expr=w23>=pi23U*d23+pi23*d23U-pi23U*d23U)
    model.ineWUG31=Constraint(expr=w31>=pi31U*d31+pi31*d31U-pi31U*d31U)
    model.ineWUG32=Constraint(expr=w32>=pi32U*d32+pi32*d32U-pi32U*d32U)
    model.ineWUG33=Constraint(expr=w33>=pi33U*d33+pi33*d33U-pi33U*d33U)

    model.ineWLL11=Constraint(expr=w11<=pi11U*d11+pi11*d11L-pi11U*d11L)
    model.ineWLL12=Constraint(expr=w12<=pi12U*d12+pi12*d12L-pi12U*d12L)
    model.ineWLL13=Constraint(expr=w13<=pi13U*d13+pi13*d13L-pi13U*d13L)
    model.ineWLL21=Constraint(expr=w21<=pi21U*d21+pi21*d21L-pi21U*d21L)
    model.ineWLL22=Constraint(expr=w22<=pi22U*d22+pi22*d22L-pi22U*d22L)
    model.ineWLL23=Constraint(expr=w23<=pi23U*d23+pi23*d23L-pi23U*d23L)
    model.ineWLL31=Constraint(expr=w31<=pi31U*d31+pi31*d31L-pi31U*d31L)
    model.ineWLL32=Constraint(expr=w32<=pi32U*d32+pi32*d32L-pi32U*d32L)
    model.ineWLL33=Constraint(expr=w33<=pi33U*d33+pi33*d33L-pi33U*d33L)

    model.ineWUL11=Constraint(expr=w11<=pi11L*d11+pi11*d11U-pi11L*d11U)
    model.ineWUL12=Constraint(expr=w12<=pi12L*d12+pi12*d12U-pi12L*d12U)
    model.ineWUL13=Constraint(expr=w13<=pi13L*d13+pi13*d13U-pi13L*d13U)
    model.ineWUL21=Constraint(expr=w21<=pi21L*d21+pi21*d21U-pi21L*d21U)
    model.ineWUL22=Constraint(expr=w22<=pi22L*d22+pi22*d22U-pi22L*d22U)
    model.ineWUL23=Constraint(expr=w23<=pi23L*d23+pi23*d23U-pi23L*d23U)
    model.ineWUL31=Constraint(expr=w31<=pi31L*d31+pi31*d31U-pi31L*d31U)
    model.ineWUL32=Constraint(expr=w32<=pi32L*d32+pi32*d32U-pi32L*d32U)
    model.ineWUL33=Constraint(expr=w33<=pi33L*d33+pi33*d33U-pi33L*d33U)

    model.distanceC = Constraint(expr = w11+w12+w13+w21+w22+w23+w31+w32+w33 <= distanceThreshold)

    # x1=x1Solution
    # x2=x2Solution
    # x1=xVectorSection[0]
    # x2=xVectorSection[1]
    # x3=xVectorSection[2]
    # x4=xVectorSection[3]
    # x5=xVectorSection[4]
    # x6=xVectorSection[5]
    # x7=xVectorSection[6]
    # x8=xVectorSection[7]
    # x9=xVectorSection[8]
    # x10=xVectorSection[9]

    xSquare=lamdaWeight*sum(sigmaWeightMatrix[i]*xVectorSection[i]**2 for i in model.indices)

    #objfunc=-(wLH1*sum( miuLH1[i]*xVectorSection[i] for i in model.indices ) + wLH2*sum( miuLH2[i]*xVectorSection[i] for i in model.indices ) + wLH3*sum( miuLH3[i]*xVectorSection[i] for i in model.indices)) + sum(xVectorSection[i]**2 for i in model.indices)
    objfunc=-(wLH1*sum( miuLH1[i]*xVectorSection[i] for i in model.indices ) + wLH2*sum( miuLH2[i]*xVectorSection[i] for i in model.indices ) + wLH3*sum( miuLH3[i]*xVectorSection[i] for i in model.indices)) + xSquare

    model.obj1 = Objective(expr=objfunc, sense=maximize)
    opt=SolverFactory('ipopt', executable='/content/ipopt')
    #opt=SolverFactory('gurobi', solver_io="python")
    results = opt.solve(model, tee=False) # solves and updates instance
    print('\\\\PrimalObj= ', model.obj1())
    print('\\\\PrimalmiuLH1',[model.miuLH1[i].value for i in model.indices])
    print('\\\\PrimalalphaLH1',[model.sigmaLH1[i].value for i in model.indices])
    print('\\\\PrimalmiuLH2',[model.miuLH2[i].value for i in model.indices])
    print('\\\\PrimalalphaLH2',[model.sigmaLH2[i].value for i in model.indices])
    print('\\\\PrimalmiuLH3',[model.miuLH3[i].value for i in model.indices])
    print('\\\\PrimalalphaLH3',[model.sigmaLH3[i].value for i in model.indices])

    # print('\\\\PrimalmiuLH1',model.miuLH1())
    # print('\\\\PrimalalphaLH1',model.sigmaLH1())
    # print('\\\\PrimalmiuLH2',model.miuLH2())
    # print('\\\\PrimalalphaLH2',model.sigmaLH2())
    # print('\\\\PrimalmiuLH3',model.miuLH3())
    # print('\\\\PrimalalphaLH3',model.sigmaLH3())

    # print(pi11())
    # print(pi12())
    # print(pi13())
    # print(pi21())
    # print(pi22())
    # print(pi23())
    # print(pi31())
    # print(pi32())
    # print(pi33())
    # pi11n=pi11()
    # pi12n=pi12()
    # pi13n=pi13()
    # pi21n=pi21()
    # pi22n=pi22()
    # pi23n=pi23()
    # pi31n=pi31()
    # pi32n=pi32()
    # pi33n=pi33()

    zeroPrevent=10E-3
    max_McGap=max(abs((w11()-pi11()*d11())/(pi11()*d11()+zeroPrevent)),abs((w12()-pi12()*d12())/(pi12()*d12()+zeroPrevent)),abs((w13()-pi13()*d13())/(pi13()*d13()+zeroPrevent)),
                  abs((w21()-pi21()*d21())/(pi21()*d21()+zeroPrevent)),abs((w22()-pi22()*d22())/(pi22()*d22()+zeroPrevent)),abs((w23()-pi23()*d23())/(pi23()*d23()+zeroPrevent)),
                  abs((w31()-pi31()*d31())/(pi31()*d31()+zeroPrevent)),abs((w32()-pi32()*d32())/(pi32()*d32()+zeroPrevent)),abs((w33()-pi33()*d33())/(pi33()*d33()+zeroPrevent)))

    # zeroPrevent=10E-8
    # max_McGap=max(abs((w11()-pi11()*d11())/(w11()+zeroPrevent)),abs((w12()-pi12()*d12())/(w12()+zeroPrevent)),abs((w13()-pi13()*d13())/(w13()+zeroPrevent)),
    #               abs((w21()-pi21()*d21())/(w21()+zeroPrevent)),abs((w22()-pi22()*d22())/(w22()+zeroPrevent)),abs((w23()-pi23()*d23())/(w23()+zeroPrevent)),
    #               abs((w31()-pi31()*d31())/(w31()+zeroPrevent)),abs((w32()-pi32()*d32())/(w32()+zeroPrevent)),abs((w33()-pi33()*d33())/(w33()+zeroPrevent)))

    # max_McGap=max(abs(w11()-pi11()*d11()),abs(w12()-pi12()*d12()),abs(w13()-pi13()*d13()),
    #               abs(w21()-pi21()*d21()),abs(w22()-pi22()*d22()),abs(w23()-pi23()*d23()),
    #               abs(w31()-pi31()*d31()),abs(w32()-pi32()*d32()),abs(w33()-pi33()*d33()))

    #print('max_McGap',max_McGap)

    # return pi11n,pi12n,pi13n,pi21n,pi22n,pi23n,pi31n,pi32n,pi33n
    return max_McGap


############################################################################### for loop
#distanceThreshold=10000
# #firstUB=distanceThreshold*1000
# firstUB=10
# dUpperBound=[firstUB,firstUB,firstUB,firstUB,firstUB,firstUB,firstUB,firstUB,firstUB]
# #dUpperBound=[0,firstUB,firstUB,firstUB,firstUB,0,firstUB,firstUB,0]

mccormickGapList=[]
objList=[]


wUH1=0.126
wUH2=0.622
wUH3=0.252
miuUH1=[4.41, 8.75, 12.98, 17.54, 20.77, 25.52, 29.71, 34.27, 38.40, 42.94]
miuUH2=[0.88, 1.43, 2.28, 3.11, 3.85, 4.42, 5.37, 5.90, 6.97, 7.25]
miuUH3=[2.60, 5.16, 7.69, 9.86, 12.91, 15.14, 17.49, 20.37, 22.33, 24.99]
muZeroVector=[miuUH1,miuUH2,miuUH3]

##############################
wLH1=0.126
wLH2=0.622
wLH3=0.252

##协方差
sigmaUH1=[1.94, 3.31, 6.67, 6.59, 8.26, 9.97, 11.40, 12.94, 13.73, 20.22]
sigmaUH2=[1.92, 3.63, 5.27, 7.26, 8.79, 12.12, 11.66, 14.56, 15.50, 15.69]
sigmaUH3=[1.93, 3.88, 6.43, 7.70, 7.65, 10.84, 13.48, 14.14, 20.26, 21.38]

#########################
sigmaWeightMatrix=[3.221418451485, 10.821902782519, 20.890576417372, 35.918669993320, 55.006095785363, 74.961495975950, 101.656052942563, 122.808065450563, 154.652601762192, 191.586515377861]
lamdaWeight=0.05
#########################

def normFunction(muZero_m,muJ):
  result = [a - b for a, b in zip(muZero_m, muJ)]
  return sum([normElement**2 for normElement in result])

def functionG_Muj(Weight_j,mu_j,x,gamma_j,muZero):
  WjmujTx=-(Weight_j*sum([mu_j_element * x_element for mu_j_element, x_element in zip(mu_j, x)]))
  gammaNorm=sum([gamma_j_element*normFunction(muZero_element,mu_j) for gamma_j_element,muZero_element in zip(gamma_j,muZero)])
  return WjmujTx-gammaNorm

def opFunc(iterationNum,boundUpperListF,boundLowerListF):
    d11U=boundUpperListF[0]
    d12U=boundUpperListF[1]
    d13U=boundUpperListF[2]
    d21U=boundUpperListF[3]
    d22U=boundUpperListF[4]
    d23U=boundUpperListF[5]
    d31U=boundUpperListF[6]
    d32U=boundUpperListF[7]
    d33U=boundUpperListF[8]

    d11L=boundLowerListF[0]
    d12L=boundLowerListF[1]
    d13L=boundLowerListF[2]
    d21L=boundLowerListF[3]
    d22L=boundLowerListF[4]
    d23L=boundLowerListF[5]
    d31L=boundLowerListF[6]
    d32L=boundLowerListF[7]
    d33L=boundLowerListF[8]

    model = ConcreteModel()

    model.lamda1=Var(within=Reals,initialize=0.1)
    model.lamda2=Var(within=Reals,initialize=0.1)
    model.lamda3=Var(within=Reals,initialize=0.1)
    model.lamda4=Var(within=Reals,initialize=0.1)
    model.lamda5=Var(within=Reals,initialize=0.1)
    model.lamda6=Var(within=Reals,initialize=0.1)

    lamda1=model.lamda1
    lamda2=model.lamda2
    lamda3=model.lamda3
    lamda4=model.lamda4
    lamda5=model.lamda5
    lamda6=model.lamda6

    model.r11=Var(within=NonNegativeReals,  initialize=0.0)
    model.r12=Var(within=NonNegativeReals,  initialize=0.0)
    model.r13=Var(within=NonNegativeReals,  initialize=0.0)
    model.r21=Var(within=NonNegativeReals,  initialize=0.0)
    model.r22=Var(within=NonNegativeReals,  initialize=0.0)
    model.r23=Var(within=NonNegativeReals,  initialize=0.0)
    model.r31=Var(within=NonNegativeReals,  initialize=0.0)
    model.r32=Var(within=NonNegativeReals,  initialize=0.0)
    model.r33=Var(within=NonNegativeReals,  initialize=0.0)

    r11=model.r11
    r12=model.r12
    r13=model.r13
    r21=model.r21
    r22=model.r22
    r23=model.r23
    r31=model.r31
    r32=model.r32
    r33=model.r33

    model.n11L=Var(within=NonNegativeReals,  initialize=0.0)
    model.n12L=Var(within=NonNegativeReals,  initialize=0.0)
    model.n13L=Var(within=NonNegativeReals,  initialize=0.0)
    model.n21L=Var(within=NonNegativeReals,  initialize=0.0)
    model.n22L=Var(within=NonNegativeReals,  initialize=0.0)
    model.n23L=Var(within=NonNegativeReals,  initialize=0.0)
    model.n31L=Var(within=NonNegativeReals,  initialize=0.0)
    model.n32L=Var(within=NonNegativeReals,  initialize=0.0)
    model.n33L=Var(within=NonNegativeReals,  initialize=0.0)

    n11L=model.n11L
    n12L=model.n12L
    n13L=model.n13L
    n21L=model.n21L
    n22L=model.n22L
    n23L=model.n23L
    n31L=model.n31L
    n32L=model.n32L
    n33L=model.n33L

    model.n11U=Var(within=NonNegativeReals,  initialize=0.0)
    model.n12U=Var(within=NonNegativeReals,  initialize=0.0)
    model.n13U=Var(within=NonNegativeReals,  initialize=0.0)
    model.n21U=Var(within=NonNegativeReals,  initialize=0.0)
    model.n22U=Var(within=NonNegativeReals,  initialize=0.0)
    model.n23U=Var(within=NonNegativeReals,  initialize=0.0)
    model.n31U=Var(within=NonNegativeReals,  initialize=0.0)
    model.n32U=Var(within=NonNegativeReals,  initialize=0.0)
    model.n33U=Var(within=NonNegativeReals,  initialize=0.0)

    n11U=model.n11U
    n12U=model.n12U
    n13U=model.n13U
    n21U=model.n21U
    n22U=model.n22U
    n23U=model.n23U
    n31U=model.n31U
    n32U=model.n32U
    n33U=model.n33U


    model.n11UHL=Var(within=NonNegativeReals,  initialize=0.0)
    model.n12UHL=Var(within=NonNegativeReals,  initialize=0.0)
    model.n13UHL=Var(within=NonNegativeReals,  initialize=0.0)
    model.n21UHL=Var(within=NonNegativeReals,  initialize=0.0)
    model.n22UHL=Var(within=NonNegativeReals,  initialize=0.0)
    model.n23UHL=Var(within=NonNegativeReals,  initialize=0.0)
    model.n31UHL=Var(within=NonNegativeReals,  initialize=0.0)
    model.n32UHL=Var(within=NonNegativeReals,  initialize=0.0)
    model.n33UHL=Var(within=NonNegativeReals,  initialize=0.0)

    n11UHL=model.n11UHL
    n12UHL=model.n12UHL
    n13UHL=model.n13UHL
    n21UHL=model.n21UHL
    n22UHL=model.n22UHL
    n23UHL=model.n23UHL
    n31UHL=model.n31UHL
    n32UHL=model.n32UHL
    n33UHL=model.n33UHL

    model.n11UHU=Var(within=NonNegativeReals,  initialize=0.0)
    model.n12UHU=Var(within=NonNegativeReals,  initialize=0.0)
    model.n13UHU=Var(within=NonNegativeReals,  initialize=0.0)
    model.n21UHU=Var(within=NonNegativeReals,  initialize=0.0)
    model.n22UHU=Var(within=NonNegativeReals,  initialize=0.0)
    model.n23UHU=Var(within=NonNegativeReals,  initialize=0.0)
    model.n31UHU=Var(within=NonNegativeReals,  initialize=0.0)
    model.n32UHU=Var(within=NonNegativeReals,  initialize=0.0)
    model.n33UHU=Var(within=NonNegativeReals,  initialize=0.0)

    n11UHU=model.n11UHU
    n12UHU=model.n12UHU
    n13UHU=model.n13UHU
    n21UHU=model.n21UHU
    n22UHU=model.n22UHU
    n23UHU=model.n23UHU
    n31UHU=model.n31UHU
    n32UHU=model.n32UHU
    n33UHU=model.n33UHU


    model.theta=Var(within=NonNegativeReals,  initialize=0.0)
    theta=model.theta


    # wUH1=0.6
    # wUH2=0.25
    # wUH3=0.15
    # miuUH1=[1*0.75,2*0.75,3*0.75,4*0.75,5*0.75,6*0.75,7*0.75,8*0.75,9*0.75,10*0.75]
    # miuUH2=[1*2.5,2*2.5,3*2.5,4*2.5,5*2.5,6*2.5,7*2.5,8*2.5,9*2.5,10*2.5]
    # miuUH3=[1*4.25,2*4.25,3*4.25,4*4.25,5*4.25,6*4.25,7*4.25,8*4.25,9*4.25,10*4.25]
    # muZeroVector=[miuUH1,miuUH2,miuUH3]

    # ##############################
    # wLH1=0.6
    # wLH2=0.25
    # wLH3=0.15
    ##############################
    model.x1=Var(within=NonNegativeReals, initialize=0.1)
    model.x2=Var(within=NonNegativeReals, initialize=0.1)
    model.x3=Var(within=NonNegativeReals, initialize=0.1)
    model.x4=Var(within=NonNegativeReals, initialize=0.1)
    model.x5=Var(within=NonNegativeReals, initialize=0.1)
    model.x6=Var(within=NonNegativeReals, initialize=0.1)
    model.x7=Var(within=NonNegativeReals, initialize=0.1)
    model.x8=Var(within=NonNegativeReals, initialize=0.1)
    model.x9=Var(within=NonNegativeReals, initialize=0.1)
    model.x10=Var(within=NonNegativeReals, initialize=0.1)
    x1=model.x1
    x2=model.x2
    x3=model.x3
    x4=model.x4
    x5=model.x5
    x6=model.x6
    x7=model.x7
    x8=model.x8
    x9=model.x9
    x10=model.x10
    xVector=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10]
    gammaVector_1=[r11,r21,r31]
    gammaVector_2=[r12,r22,r32]
    gammaVector_3=[r13,r23,r33]

    #obj_miu1Star=wLH1/(2*r11+2*r21+2*r31)*xVector+(r11*miuUH1+r21*miuUH2+r31*miuUH3)/(r11+r21+r31)    [r11 * value for value in miuUH1]
    #obj_miu1Star=[-wLH1/(2*r11+2*r21+2*r31)*xVectorElement for xVectorElement in xVector]+(r11*np.array(miuUH1)+r21*np.array(miuUH2)+r31*np.array(miuUH3))/(r11+r21+r31)
    obj_miu1Star=[-wLH1/(2*r11+2*r21+2*r31)*xVectorElement for xVectorElement in xVector]+\
        (np.array([r11, r11*2, r11*3])+[r21, r21*2, r21*3]+[r31, r31*2, r31*3])/(r11+r21+r31)
    obj_gmiu1Star=functionG_Muj(wLH1,obj_miu1Star,xVector,gammaVector_1,muZeroVector)

    #obj_miu2Star=wLH2/(2*r12+2*r22+2*r32)*xVector+(r12*miuUH1+r22*miuUH2+r32*miuUH3)/(r12+r22+r32)
    obj_miu2Star=[-wLH2/(2*r12+2*r22+2*r32)*xVectorElement for xVectorElement in xVector]+(r12*np.array(miuUH1)+r22*np.array(miuUH2)+r32*np.array(miuUH3))/(r12+r22+r32)
    obj_gmiu2Star=functionG_Muj(wLH2,obj_miu2Star,xVector,gammaVector_2,muZeroVector)

    #obj_miu3Star=wLH3/(2*r13+2*r23+2*r33)*xVector+(r13*miuUH1+r23*miuUH2+r33*miuUH3)/(r13+r23+r33)
    obj_miu3Star=[-wLH3/(2*r13+2*r23+2*r33)*xVectorElement for xVectorElement in xVector]+(r13*np.array(miuUH1)+r23*np.array(miuUH2)+r33*np.array(miuUH3))/(r13+r23+r33)
    obj_gmiu3Star=functionG_Muj(wLH3,obj_miu3Star,xVector,gammaVector_3,muZeroVector)

    # ##协方差
    # sigmaUH1=[1*1.8,2*1.8,3*1.8,4*1.8,5*1.8,6*1.8,7*1.8,8*1.8,9*1.8,10*1.8]
    # sigmaUH2=[1*1.8,2*1.8,3*1.8,4*1.8,5*1.8,6*1.8,7*1.8,8*1.8,9*1.8,10*1.8]
    # sigmaUH3=[1*1.8,2*1.8,3*1.8,4*1.8,5*1.8,6*1.8,7*1.8,8*1.8,9*1.8,10*1.8]

    sigma1DiagHalf=[r11*sigmaUH1Element**0.5+r21*sigmaUH2Element**0.5+r31*sigmaUH3Element**0.5 for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    traceSigma1Square=sum([sigma1DiagHalfElement**2 for sigma1DiagHalfElement in sigma1DiagHalf])
    sigma1Diag=[r11*sigmaUH1Element+r21*sigmaUH2Element+r31*sigmaUH3Element for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    traceSigma1minus=sum([-sigma1DiagElement for sigma1DiagElement in sigma1Diag])
    obj_sigma1Star=(1/(r11+r21+r31))*traceSigma1Square+traceSigma1minus


    sigma2DiagHalf=[r12*sigmaUH1Element**0.5+r22*sigmaUH2Element**0.5+r32*sigmaUH3Element**0.5 for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    traceSigma2Square=sum([sigma2DiagHalfElement**2 for sigma2DiagHalfElement in sigma2DiagHalf])
    sigma2Diag=[r12*sigmaUH1Element+r22*sigmaUH2Element+r32*sigmaUH3Element for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    traceSigma2minus=sum([-sigma2DiagElement for sigma2DiagElement in sigma2Diag])
    obj_sigma2Star=(1/(r12+r22+r32))*traceSigma2Square+traceSigma2minus

    sigma3DiagHalf=[r13*sigmaUH1Element**0.5+r23*sigmaUH2Element**0.5+r33*sigmaUH3Element**0.5 for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    traceSigma3Square=sum([sigma3DiagHalfElement**2 for sigma3DiagHalfElement in sigma3DiagHalf])
    sigma3Diag=[r13*sigmaUH1Element+r23*sigmaUH2Element+r33*sigmaUH3Element for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    traceSigma3minus=sum([-sigma3DiagElement for sigma3DiagElement in sigma3Diag])
    obj_sigma3Star=(1/(r13+r23+r33))*traceSigma3Square+traceSigma3minus

##############################

    pi11L=resultsRange['pi11'][0]
    pi12L=resultsRange['pi12'][0]
    pi13L=resultsRange['pi13'][0]
    pi21L=resultsRange['pi21'][0]
    pi22L=resultsRange['pi22'][0]
    pi23L=resultsRange['pi23'][0]
    pi31L=resultsRange['pi31'][0]
    pi32L=resultsRange['pi32'][0]
    pi33L=resultsRange['pi33'][0]

    pi11U=resultsRange['pi11'][1]
    pi12U=resultsRange['pi12'][1]
    pi13U=resultsRange['pi13'][1]
    pi21U=resultsRange['pi21'][1]
    pi22U=resultsRange['pi22'][1]
    pi23U=resultsRange['pi23'][1]
    pi31U=resultsRange['pi31'][1]
    pi32U=resultsRange['pi32'][1]
    pi33U=resultsRange['pi33'][1]

    ############## obj ##################
    obj_lamdaTW=wUH1*lamda1+wUH2*lamda2+wUH3*lamda3+wLH1*lamda4+wLH2*lamda5+wLH3*lamda6
    obj_lamdathetaDistance=theta*distanceThreshold
    obj_nLTpiLdL=n11L*pi11L*d11L+n12L*pi12L*d12L+n13L*pi13L*d13L+n21L*pi21L*d21L+n22L*pi22L*d22L+n23L*pi23L*d23L+n31L*pi31L*d31L+n32L*pi32L*d32L+n33L*pi33L*d33L
    obj_nUTpiUdU=n11U*pi11U*d11U+n12U*pi12U*d12U+n13U*pi13U*d13U+n21U*pi21U*d21U+n22U*pi22U*d22U+n23U*pi23U*d23U+n31U*pi31U*d31U+n32U*pi32U*d32U+n33U*pi33U*d33U

    obj_nUHLTpiLdL=-n11UHL*pi11L*d11U-n12UHL*pi12L*d12U-n13UHL*pi13L*d13U-n21UHL*pi21L*d21U-n22UHL*pi22L*d22U-n23UHL*pi23L*d23U-n31UHL*pi31L*d31U-n32UHL*pi32L*d32U-n33UHL*pi33L*d33U
    obj_nUHUTpiUdU=-n11UHU*pi11U*d11L-n12UHU*pi12U*d12L-n13UHU*pi13U*d13L-n21UHU*pi21U*d21L-n22UHU*pi22U*d22L-n23UHU*pi23U*d23L-n31UHU*pi31U*d31L-n32UHU*pi32U*d32L-n33UHU*pi33U*d33L

    # obj_rTalpha=-(r11*alphaUH1**2+r12*alphaUH1**2+r13*alphaUH1**2+r21*alphaUH2**2+r22*alphaUH2**2+r23*alphaUH2**2+r31*alphaUH3**2+r32*alphaUH3**2+r33*alphaUH3**2)

    # obj_miu1Star=( (wLH1**2*x1**2)/4 + (miuUH1*r11+miuUH2*r21+miuUH3*r31)*wLH1*x1 - ( (r11*r21*(miuUH1-miuUH2)**2) + (r11*r31*(miuUH1-miuUH3)**2) + (r21*r31*(miuUH2-miuUH3)**2) ) )/(r11+r21+r31)
    # obj_miu2Star=( (wLH2**2*x1**2)/4 + (miuUH1*r12+miuUH2*r22+miuUH3*r32)*wLH2*x1 - ( (r12*r22*(miuUH1-miuUH2)**2) + (r12*r32*(miuUH1-miuUH3)**2) + (r22*r32*(miuUH2-miuUH3)**2) ) )/(r12+r22+r32)
    # obj_miu3Star=( (wLH3**2*x1**2)/4 + (miuUH1*r13+miuUH2*r23+miuUH3*r33)*wLH3*x1 - ( (r13*r23*(miuUH1-miuUH2)**2) + (r13*r33*(miuUH1-miuUH3)**2) + (r23*r33*(miuUH2-miuUH3)**2) ) )/(r13+r23+r33)
    # obj_alpha1Star=(r11*alphaUH1+r21*alphaUH2+r31*alphaUH3)**2/(r11+r21+r31)
    # obj_alpha2Star=(r12*alphaUH1+r22*alphaUH2+r32*alphaUH3)**2/(r12+r22+r32)
    # obj_alpha3Star=(r13*alphaUH1+r23*alphaUH2+r33*alphaUH3)**2/(r13+r23+r33)

    functoinGandfunctionH=obj_gmiu1Star+obj_gmiu2Star+obj_gmiu3Star+obj_sigma1Star+obj_sigma2Star+obj_sigma3Star


    xSquare=lamdaWeight*sum([sigmaWeight*xElement**2 for xElement,sigmaWeight in zip(xVector,sigmaWeightMatrix)])

    model.obj1 = Objective(expr=obj_lamdaTW+obj_lamdathetaDistance+obj_nLTpiLdL+obj_nUTpiUdU+obj_nUHLTpiLdL+obj_nUHUTpiUdU+functoinGandfunctionH+xSquare, sense=minimize)
############ constraint ##############

    #########C1#############

    C1_T1_1=-(lamda1+lamda4)
    C1_T1_2=-(lamda1+lamda5)
    C1_T1_3=-(lamda1+lamda6)
    C1_T1_4=-(lamda2+lamda4)
    C1_T1_5=-(lamda2+lamda5)
    C1_T1_6=-(lamda2+lamda6)
    C1_T1_7=-(lamda3+lamda4)
    C1_T1_8=-(lamda3+lamda5)
    C1_T1_9=-(lamda3+lamda6)


    C1_T2_1=-(n11L*d11L)
    C1_T2_2=-(n12L*d12L)
    C1_T2_3=-(n13L*d13L)
    C1_T2_4=-(n21L*d21L)
    C1_T2_5=-(n22L*d22L)
    C1_T2_6=-(n23L*d23L)
    C1_T2_7=-(n31L*d31L)
    C1_T2_8=-(n32L*d32L)
    C1_T2_9=-(n33L*d33L)

    C1_T3_1=-(n11U*d11U)
    C1_T3_2=-(n12U*d12U)
    C1_T3_3=-(n13U*d13U)
    C1_T3_4=-(n21U*d21U)
    C1_T3_5=-(n22U*d22U)
    C1_T3_6=-(n23U*d23U)
    C1_T3_7=-(n31U*d31U)
    C1_T3_8=-(n32U*d32U)
    C1_T3_9=-(n33U*d33U)

    C1_T4_1=(n11UHL*d11U)
    C1_T4_2=(n12UHL*d12U)
    C1_T4_3=(n13UHL*d13U)
    C1_T4_4=(n21UHL*d21U)
    C1_T4_5=(n22UHL*d22U)
    C1_T4_6=(n23UHL*d23U)
    C1_T4_7=(n31UHL*d31U)
    C1_T4_8=(n32UHL*d32U)
    C1_T4_9=(n33UHL*d33U)

    C1_T5_1=(n11UHU*d11L)
    C1_T5_2=(n12UHU*d12L)
    C1_T5_3=(n13UHU*d13L)
    C1_T5_4=(n21UHU*d21L)
    C1_T5_5=(n22UHU*d22L)
    C1_T5_6=(n23UHU*d23L)
    C1_T5_7=(n31UHU*d31L)
    C1_T5_8=(n32UHU*d32L)
    C1_T5_9=(n33UHU*d33L)

    model.C1_1 = Constraint(expr = C1_T1_1+C1_T2_1+C1_T3_1+C1_T4_1+C1_T5_1<=0.0)
    model.C1_2 = Constraint(expr = C1_T1_2+C1_T2_2+C1_T3_2+C1_T4_2+C1_T5_2<=0.0)
    model.C1_3 = Constraint(expr = C1_T1_3+C1_T2_3+C1_T3_3+C1_T4_3+C1_T5_3<=0.0)
    model.C1_4 = Constraint(expr = C1_T1_4+C1_T2_4+C1_T3_4+C1_T4_4+C1_T5_4<=0.0)
    model.C1_5 = Constraint(expr = C1_T1_5+C1_T2_5+C1_T3_5+C1_T4_5+C1_T5_5<=0.0)
    model.C1_6 = Constraint(expr = C1_T1_6+C1_T2_6+C1_T3_6+C1_T4_6+C1_T5_6<=0.0)
    model.C1_7 = Constraint(expr = C1_T1_7+C1_T2_7+C1_T3_7+C1_T4_7+C1_T5_7<=0.0)
    model.C1_8 = Constraint(expr = C1_T1_8+C1_T2_8+C1_T3_8+C1_T4_8+C1_T5_8<=0.0)
    model.C1_9 = Constraint(expr = C1_T1_9+C1_T2_9+C1_T3_9+C1_T4_9+C1_T5_9<=0.0)

    ########C2###############
    C2_T1_1=-theta
    C2_T1_2=-theta
    C2_T1_3=-theta
    C2_T1_4=-theta
    C2_T1_5=-theta
    C2_T1_6=-theta
    C2_T1_7=-theta
    C2_T1_8=-theta
    C2_T1_9=-theta

    C2_T2_1=n11L
    C2_T2_2=n12L
    C2_T2_3=n13L
    C2_T2_4=n21L
    C2_T2_5=n22L
    C2_T2_6=n23L
    C2_T2_7=n31L
    C2_T2_8=n32L
    C2_T2_9=n33L

    C2_T3_1=n11U
    C2_T3_2=n12U
    C2_T3_3=n13U
    C2_T3_4=n21U
    C2_T3_5=n22U
    C2_T3_6=n23U
    C2_T3_7=n31U
    C2_T3_8=n32U
    C2_T3_9=n33U

    C2_T4_1=-n11UHL
    C2_T4_2=-n12UHL
    C2_T4_3=-n13UHL
    C2_T4_4=-n21UHL
    C2_T4_5=-n22UHL
    C2_T4_6=-n23UHL
    C2_T4_7=-n31UHL
    C2_T4_8=-n32UHL
    C2_T4_9=-n33UHL

    C2_T5_1=-n11UHU
    C2_T5_2=-n12UHU
    C2_T5_3=-n13UHU
    C2_T5_4=-n21UHU
    C2_T5_5=-n22UHU
    C2_T5_6=-n23UHU
    C2_T5_7=-n31UHU
    C2_T5_8=-n32UHU
    C2_T5_9=-n33UHU

    model.C2_1 = Constraint(expr = C2_T1_1+C2_T2_1+C2_T3_1+C2_T4_1+C2_T5_1<=0.0)
    model.C2_2 = Constraint(expr = C2_T1_2+C2_T2_2+C2_T3_2+C2_T4_2+C2_T5_2<=0.0)
    model.C2_3 = Constraint(expr = C2_T1_3+C2_T2_3+C2_T3_3+C2_T4_3+C2_T5_3<=0.0)
    model.C2_4 = Constraint(expr = C2_T1_4+C2_T2_4+C2_T3_4+C2_T4_4+C2_T5_4<=0.0)
    model.C2_5 = Constraint(expr = C2_T1_5+C2_T2_5+C2_T3_5+C2_T4_5+C2_T5_5<=0.0)
    model.C2_6 = Constraint(expr = C2_T1_6+C2_T2_6+C2_T3_6+C2_T4_6+C2_T5_6<=0.0)
    model.C2_7 = Constraint(expr = C2_T1_7+C2_T2_7+C2_T3_7+C2_T4_7+C2_T5_7<=0.0)
    model.C2_8 = Constraint(expr = C2_T1_8+C2_T2_8+C2_T3_8+C2_T4_8+C2_T5_8<=0.0)
    model.C2_9 = Constraint(expr = C2_T1_9+C2_T2_9+C2_T3_9+C2_T4_9+C2_T5_9<=0.0)

    ########C3###############
    C3_T1_1=r11
    C3_T1_2=r12
    C3_T1_3=r13
    C3_T1_4=r21
    C3_T1_5=r22
    C3_T1_6=r23
    C3_T1_7=r31
    C3_T1_8=r32
    C3_T1_9=r33

    C3_T2_1=-(n11L*pi11L)
    C3_T2_2=-(n12L*pi12L)
    C3_T2_3=-(n13L*pi13L)
    C3_T2_4=-(n21L*pi12L)
    C3_T2_5=-(n22L*pi22L)
    C3_T2_6=-(n23L*pi23L)
    C3_T2_7=-(n31L*pi31L)
    C3_T2_8=-(n32L*pi32L)
    C3_T2_9=-(n33L*pi33L)

    C3_T3_1=-(n11U*pi11U)
    C3_T3_2=-(n12U*pi12U)
    C3_T3_3=-(n13U*pi13U)
    C3_T3_4=-(n21U*pi12U)
    C3_T3_5=-(n22U*pi22U)
    C3_T3_6=-(n23U*pi23U)
    C3_T3_7=-(n31U*pi31U)
    C3_T3_8=-(n32U*pi32U)
    C3_T3_9=-(n33U*pi33U)

    C3_T4_1=(n11UHL*pi11L)
    C3_T4_2=(n12UHL*pi12L)
    C3_T4_3=(n13UHL*pi13L)
    C3_T4_4=(n21UHL*pi12L)
    C3_T4_5=(n22UHL*pi22L)
    C3_T4_6=(n23UHL*pi23L)
    C3_T4_7=(n31UHL*pi31L)
    C3_T4_8=(n32UHL*pi32L)
    C3_T4_9=(n33UHL*pi33L)

    C3_T5_1=(n11UHU*pi11U)
    C3_T5_2=(n12UHU*pi12U)
    C3_T5_3=(n13UHU*pi13U)
    C3_T5_4=(n21UHU*pi12U)
    C3_T5_5=(n22UHU*pi22U)
    C3_T5_6=(n23UHU*pi23U)
    C3_T5_7=(n31UHU*pi31U)
    C3_T5_8=(n32UHU*pi32U)
    C3_T5_9=(n33UHU*pi33U)

    model.C3_1 = Constraint(expr = C3_T1_1+C3_T2_1+C3_T3_1+C3_T4_1+C3_T5_1<=0.0)
    model.C3_2 = Constraint(expr = C3_T1_2+C3_T2_2+C3_T3_2+C3_T4_2+C3_T5_2<=0.0)
    model.C3_3 = Constraint(expr = C3_T1_3+C3_T2_3+C3_T3_3+C3_T4_3+C3_T5_3<=0.0)
    model.C3_4 = Constraint(expr = C3_T1_4+C3_T2_4+C3_T3_4+C3_T4_4+C3_T5_4<=0.0)
    model.C3_5 = Constraint(expr = C3_T1_5+C3_T2_5+C3_T3_5+C3_T4_5+C3_T5_5<=0.0)
    model.C3_6 = Constraint(expr = C3_T1_6+C3_T2_6+C3_T3_6+C3_T4_6+C3_T5_6<=0.0)
    model.C3_7 = Constraint(expr = C3_T1_7+C3_T2_7+C3_T3_7+C3_T4_7+C3_T5_7<=0.0)
    model.C3_8 = Constraint(expr = C3_T1_8+C3_T2_8+C3_T3_8+C3_T4_8+C3_T5_8<=0.0)
    model.C3_9 = Constraint(expr = C3_T1_9+C3_T2_9+C3_T3_9+C3_T4_9+C3_T5_9<=0.0)

    ########C4###############
    model.C4=Constraint(expr = x1+x2+x3+x4+x5+x6+x7+x8+x9+x10==1)
    ###########################################################
    ###########################################################
    opt=SolverFactory('ipopt', executable='/content/ipopt')
    #opt=SolverFactory('gurobi', solver_io='python')

    results = opt.solve(model, tee=True) # solves and updates instance
    objResult=model.obj1()
    print('\nDual Profit = ', objResult)
    objList.append(objResult)
    r11Solution=model.r11()
    r12Solution=model.r12()
    r13Solution=model.r13()
    r21Solution=model.r21()
    r22Solution=model.r22()
    r23Solution=model.r23()
    r31Solution=model.r31()
    r32Solution=model.r32()
    r33Solution=model.r33()
    x1Solution=model.x1()
    x2Solution=model.x2()
    x3Solution=model.x3()
    x4Solution=model.x4()
    x5Solution=model.x5()
    x6Solution=model.x6()
    x7Solution=model.x7()
    x8Solution=model.x8()
    x9Solution=model.x9()
    x10Solution=model.x10()
    xVectorSection=[x1Solution,x2Solution,x3Solution,x4Solution,x5Solution,x6Solution,x7Solution,x8Solution,x9Solution,x10Solution]
    print("x solution",xVectorSection)

    # miuLH1=wLH1*x1Solution/(2*(r11Solution+r21Solution+r31Solution)) + (r11Solution*miuUH1+r21Solution*miuUH2+r31Solution*miuUH3)/(r11Solution+r21Solution+r31Solution)
    # miuLH2=wLH2*x1Solution/(2*(r12Solution+r22Solution+r32Solution)) + (r12Solution*miuUH1+r22Solution*miuUH2+r32Solution*miuUH3)/(r12Solution+r22Solution+r32Solution)
    # miuLH3=wLH3*x1Solution/(2*(r13Solution+r23Solution+r33Solution)) + (r13Solution*miuUH1+r23Solution*miuUH2+r33Solution*miuUH3)/(r13Solution+r23Solution+r33Solution)
    miuLH1=[-wLH1/(2*r11Solution+2*r21Solution+2*r31Solution)*xVectorElement for xVectorElement in xVectorSection]+(r11Solution*np.array(miuUH1)+r21Solution*np.array(miuUH2)+r31Solution*np.array(miuUH3))/(r11Solution+r21Solution+r31Solution)
    miuLH2=[-wLH2/(2*r12Solution+2*r22Solution+2*r32Solution)*xVectorElement for xVectorElement in xVectorSection]+(r12Solution*np.array(miuUH1)+r22Solution*np.array(miuUH2)+r32Solution*np.array(miuUH3))/(r12Solution+r22Solution+r32Solution)
    miuLH3=[-wLH3/(2*r13Solution+2*r23Solution+2*r33Solution)*xVectorElement for xVectorElement in xVectorSection]+(r13Solution*np.array(miuUH1)+r23Solution*np.array(miuUH2)+r33Solution*np.array(miuUH3))/(r13Solution+r23Solution+r33Solution)
    print('miu1',miuLH1)
    print('miu2',miuLH2)
    print('miu3',miuLH3)
    # print((wLH1*miuLH1+wLH2*miuLH2+wLH3*miuLH3)*x1Solution+x1Solution+x2Solution)

    # alphaLH1=(r11Solution*alphaUH1+r21Solution*alphaUH2+r31Solution*alphaUH3)/(r11Solution+r21Solution+r31Solution)
    # alphaLH2=(r12Solution*alphaUH1+r22Solution*alphaUH2+r32Solution*alphaUH3)/(r12Solution+r22Solution+r32Solution)
    # alphaLH3=(r13Solution*alphaUH1+r23Solution*alphaUH2+r33Solution*alphaUH3)/(r13Solution+r23Solution+r33Solution)
    sigma1DiagHalfSolution=[(r11Solution*sigmaUH1Element**0.5+r21Solution*sigmaUH2Element**0.5+r31Solution*sigmaUH3Element**0.5)/(r11Solution+r21Solution+r31Solution) for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    alphaLH1=[sigma1DiagHalfElement**2 for sigma1DiagHalfElement in sigma1DiagHalfSolution]
    sigma2DiagHalfSolution=[(r12Solution*sigmaUH1Element**0.5+r22Solution*sigmaUH2Element**0.5+r32Solution*sigmaUH3Element**0.5)/(r12Solution+r22Solution+r32Solution) for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    alphaLH2=[sigma2DiagHalfElement**2 for sigma2DiagHalfElement in sigma2DiagHalfSolution]
    sigma3DiagHalfSolution=[(r13Solution*sigmaUH1Element**0.5+r23Solution*sigmaUH2Element**0.5+r33Solution*sigmaUH3Element**0.5)/(r13Solution+r23Solution+r33Solution) for sigmaUH1Element,sigmaUH2Element,sigmaUH3Element in zip(sigmaUH1,sigmaUH2,sigmaUH3)]
    alphaLH3=[sigma3DiagHalfElement**2 for sigma3DiagHalfElement in sigma3DiagHalfSolution]
    print('alphaLH1',alphaLH1)
    print('alphaLH2',alphaLH2)
    print('alphaLH3',alphaLH3)
    ##################### primal ###################
    mccormickGap=primalFunc(xVectorSection,resultsRange,boundUpperListF,boundLowerListF)
    print('mccormickGap',mccormickGap)
    mccormickGapList.append(mccormickGap)
    ################################################改到此处；改到此处；改到此处；改到此处；改到此处；改到此处；改到此处；改到此处
    # wD11=(miuUH1-miuLH1)**2+(alphaUH1**2+alphaLH1**2-2*(alphaUH1*alphaLH1))
    # wD12=(miuUH1-miuLH2)**2+(alphaUH1**2+alphaLH2**2-2*(alphaUH1*alphaLH2))
    # wD13=(miuUH1-miuLH3)**2+(alphaUH1**2+alphaLH3**2-2*(alphaUH1*alphaLH3))
    # wD21=(miuUH2-miuLH1)**2+(alphaUH2**2+alphaLH1**2-2*(alphaUH2*alphaLH1))
    # wD22=(miuUH2-miuLH2)**2+(alphaUH2**2+alphaLH2**2-2*(alphaUH2*alphaLH2))
    # wD23=(miuUH2-miuLH3)**2+(alphaUH2**2+alphaLH3**2-2*(alphaUH2*alphaLH3))
    # wD31=(miuUH3-miuLH1)**2+(alphaUH3**2+alphaLH1**2-2*(alphaUH3*alphaLH1))
    # wD32=(miuUH3-miuLH2)**2+(alphaUH3**2+alphaLH2**2-2*(alphaUH3*alphaLH2))
    # wD33=(miuUH3-miuLH3)**2+(alphaUH3**2+alphaLH3**2-2*(alphaUH3*alphaLH3))
    wD11=normFunction(miuUH1,miuLH1) + sum([Sigma1_0Element+Sigma1Element-2*(Sigma1_0Element*Sigma1Element)**0.5 for Sigma1_0Element,Sigma1Element in zip(sigmaUH1,alphaLH1) ])
    wD12=normFunction(miuUH1,miuLH2) + sum([Sigma1_0Element+Sigma2Element-2*(Sigma1_0Element*Sigma2Element)**0.5 for Sigma1_0Element,Sigma2Element in zip(sigmaUH1,alphaLH2) ])
    wD13=normFunction(miuUH1,miuLH3) + sum([Sigma1_0Element+Sigma3Element-2*(Sigma1_0Element*Sigma3Element)**0.5 for Sigma1_0Element,Sigma3Element in zip(sigmaUH1,alphaLH3) ])
    wD21=normFunction(miuUH2,miuLH1) + sum([Sigma2_0Element+Sigma1Element-2*(Sigma2_0Element*Sigma1Element)**0.5 for Sigma2_0Element,Sigma1Element in zip(sigmaUH2,alphaLH1) ])
    wD22=normFunction(miuUH2,miuLH2) + sum([Sigma2_0Element+Sigma2Element-2*(Sigma2_0Element*Sigma2Element)**0.5 for Sigma2_0Element,Sigma2Element in zip(sigmaUH2,alphaLH2) ])
    wD23=normFunction(miuUH2,miuLH3) + sum([Sigma2_0Element+Sigma3Element-2*(Sigma2_0Element*Sigma3Element)**0.5 for Sigma2_0Element,Sigma3Element in zip(sigmaUH2,alphaLH3) ])
    wD31=normFunction(miuUH3,miuLH1) + sum([Sigma3_0Element+Sigma1Element-2*(Sigma3_0Element*Sigma1Element)**0.5 for Sigma3_0Element,Sigma1Element in zip(sigmaUH3,alphaLH1) ])
    wD32=normFunction(miuUH3,miuLH2) + sum([Sigma3_0Element+Sigma2Element-2*(Sigma3_0Element*Sigma2Element)**0.5 for Sigma3_0Element,Sigma2Element in zip(sigmaUH3,alphaLH2) ])
    wD33=normFunction(miuUH3,miuLH3) + sum([Sigma3_0Element+Sigma3Element-2*(Sigma3_0Element*Sigma3Element)**0.5 for Sigma3_0Element,Sigma3Element in zip(sigmaUH3,alphaLH3) ])
    # print(wD11)
    # print(wD12)
    # print(wD13)
    # print(wD21)
    # print(wD22)
    # print(wD23)
    # print(wD31)
    # print(wD32)
    # print(wD33)

    # parameterIncrease=1.01
    # parameterdecrease=0.99
    # parameterIncrease=1.05
    # parameterdecrease=0.95
    # parameterIncrease=2.0
    # parameterdecrease=0.0

    convergencyNum=(1000-iterationNum)/1000
    parameterIncrease=1+convergencyNum
    parameterdecrease=1-convergencyNum
    print('convergencyNum',convergencyNum)
# ##################################
    # boundIU=[wD11*parameterIncrease,wD12*parameterIncrease,wD13*parameterIncrease,wD21*parameterIncrease,wD22*parameterIncrease,wD23*parameterIncrease,wD31*parameterIncrease,wD32*parameterIncrease,wD33*parameterIncrease]
    # boundIL=[wD11*parameterdecrease,wD12*parameterdecrease,wD13*parameterdecrease,wD21*parameterdecrease,wD22*parameterdecrease,wD23*parameterdecrease,wD31*parameterdecrease,wD32*parameterdecrease,wD33*parameterdecrease]
    # print(boundIU)
    # print(boundIL)
    # return boundIU,boundIL,wD11,wD12,wD13,wD21,wD22,wD23,wD31,wD32,wD33
###################################
    boundIU=[wD11*parameterIncrease,wD12*parameterIncrease,wD13*parameterIncrease,wD21*parameterIncrease,wD22*parameterIncrease,wD23*parameterIncrease,wD31*parameterIncrease,wD32*parameterIncrease,wD33*parameterIncrease]
    boundIL=[wD11*parameterdecrease,wD12*parameterdecrease,wD13*parameterdecrease,wD21*parameterdecrease,wD22*parameterdecrease,wD23*parameterdecrease,wD31*parameterdecrease,wD32*parameterdecrease,wD33*parameterdecrease]
    c = [min(x, y) for x, y in zip(boundIU, boundUpperList_copy)]
    d = [max(x, y) for x, y in zip(boundIL, boundLowerList_copy)]
    boundIU=c
    boundIL=d
    print(boundIU)
    print(boundIL)
    return boundIU,boundIL,wD11,wD12,wD13,wD21,wD22,wD23,wD31,wD32,wD33,mccormickGap

for i in range(1000):
  [boundUpperList,boundLowerList,wD11,wD12,wD13,wD21,wD22,wD23,wD31,wD32,wD33,mccormickGapFlag]=opFunc(i,boundUpperList,boundLowerList)
  if mccormickGapFlag<0.01:
    break
