###################################################################################################
# Integrated Biological and Economic Model
# Expanded from Robert et al. (2014). Evolutionary Applications, 7:1238–1251.
# Purpose: to evaluate the benefits and potential implications of
# simultaneously using sex-specific GPM and Bt-crops
###################################################################################################
###################################
# Clear Global Environment Before Script Run
rm(list = ls())

#################################
# Load library deSolve to calculate differential equations
library(deSolve)
# Load library plot3D to create contour plots
library(plot3D)

###################################################################################################
# Create Genotype Probability Matrix
###################################################################################################
#################################
# P(i|m,n) = probability an offspring of mating between 
# type m female and type n male will be of type i

# Create empty array
K<-array(0, dim=c(3,3,3))

# 1 : KK (Homozygous Dominant)
# 2 : Kk (Heterozygous)
# 3 : kk (Homozygous Recessive)

# KK x KK -> all offspring KK
K[1,1,1]=1
K[2,1,1]=0
K[3,1,1]=0

# KK x Kk -> 0.5 Kk, 0.5 KK
K[1,1,2]=0.5
K[2,1,2]=0.5
K[3,1,2]=0

# KK X kk -> all offspring Kk
K[1,1,3]=0
K[2,1,3]=1
K[3,1,3]=0

# Kk x KK -> 0.5 Kk, 0.5 KK
K[1,2,1]=0.5
K[2,2,1]=0.5
K[3,2,1]=0

# Kk x Kk -> 0.25 KK, 0.5 Kk, 0.25 kk
K[1,2,2]=0.25
K[2,2,2]=0.5
K[3,2,2]=0.25

# Kk x kk -> 0.5 Kk, 0.5 kk
K[1,2,3]=0
K[2,2,3]=0.5
K[3,2,3]=0.5

# kk x KK -> all offspring Kk
K[1,3,1]=0
K[2,3,1]=1
K[3,3,1]=0

# kk X Kk -> 0.5 Kk, 0.5 kk
K[1,3,2]=0
K[2,3,2]=0.5
K[3,3,2]=0.5

# kk X kk -> all offspring kk
K[1,3,3]=0
K[2,3,3]=0
K[3,3,3]=1

# Apply kronecker function to acheive the P-matrix for two-locus, bi-allelic system
# Creates 9X9X9 matrix of genotype probabilities based on the following genotypes
# Note that the genotype numbers 1-9 should be consistent throughout this model
# K = dominant tg lethal allele, k = null allele/placeholder
# S = dominant Bt-susceptible allele, s = recessive Bt-resistance allele
# 1 : KKSS
# 2 : KKSs
# 3 : KKss
# 4 : KkSS
# 5 : KkSs
# 6 : Kkss
# 7 : kkSS
# 8 : kkSs
# 9 : kkss
p.matrix <- kronecker(K, K)

# # check p.matrix : if we fix first and second index (parental genotypes), then
# # summing over third index (offspring genotypes) should give one
# 
# for (j in 1:9){
#     for (k in 1:9) {
#         sum_vals=0;
#         for (i in 1:9) {
#             sum_vals=sum_vals+p.matrix[i,j,k]
#         }
#         write(sum_vals,"",append=TRUE)
#     }
# }

####################################################
#Bio parameters outside of loops

cost.s = 0.417           # cost.s = fitness cost of resistance allele (little s)
cost.K = 0.523       # cost.K = fitness cost of lethal tg allele (big K)
.leaky = 0.01           # Leakyness of tg construct on scale of 0-1, where 0 = completly effective
mort.bt = 1          # mort.bt = Mortality of juveniles on Bt - cite: Yi et al. (2015)

# Parameters set via literature
.larval.prod = 12.9     # larval.prod = average rate of larval production by females per day, lambda in paper
.dens = 2e-04           # dens = density dependence parameter, alpha in paper
.strength.dens = 3.4    # strength.dens = strength of density dependence, beta in paper
.mort.j = 0.0029        # mort.j = daily mortality rates for juveniles, mu in paper - cite: Marchioro et al (2014)
.mort.f = 0.037          # mort.f = daily mortality rates for females, mu in paper
.mort.m = 0.019          # mort.m = daily mortality rates for males, mu in paper
.emergence = 0.055       # emergence = rate of emergence to adulthood per day, v in paper
dom.K=0.505             # dom.K = dominance of the lethal tg allele
dom.s = 0               # dom.s = dominance of Bt resistance allele (recessive so 0 is realistic value)
dom.th = 1              # dom.th = dominance of Bt susceptible allele (S is dominant so 1 is default)

# Initialization of population sizes
### Different eqns for different sexes. See supp mats from Robert et al 2013. PLos One 8(9):e73233-9
equil.j <- (1/.dens)*((.mort.j+.emergence)*(((.emergence*.larval.prod)/
                                               (2*.mort.f*(.mort.j+.emergence)))-1)^(1/(.strength.dens-1)))
equil.f <- .emergence/(2*.mort.f) * equil.j
equil.m <- .emergence/(2*.mort.m) * equil.j


###############################
.release2 = 0  
.release3 = 0  
.release4 = 0 
.release5 = 0
.release6 = 0
.release7 = 0           # Change this value to release wild-type susceptible individuals
.release8 = 0
.release9 = 0

# viability = female viability coefficients, gamma in paper
# for genotypes 1 - 9
.viability1 = 0
.viability2 = 0
.viability3 = 0
.viability4 = 0
.viability5 = 0
.viability6 = 0
.viability7 = 1
.viability8 = 1
.viability9 = 1

# fitness, omega (w) in paper
# for genotypes 1 - 9
.fitness1 = (1-cost.K)
.fitness2 = (1-dom.s*cost.s)*(1-cost.K)
.fitness3 = (1-cost.s)*(1-cost.K)
.fitness4 = (1-dom.K*cost.K)    
.fitness5 = (1-dom.s*cost.s)*(1-dom.K*cost.K)
.fitness6 = (1-cost.s)*(1-dom.K*cost.K)
.fitness7 = 1
.fitness8 = (1-dom.s*cost.s)
.fitness9 = (1-cost.s)

###################################################################################################
# Loop for the NPV surface, with Refuge (q) against Releases
###################################################################################################
npv.cum.list<-array(0, dim=c(11,11,9))

init.resist.list<-c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 0.9)
#Loop over bT AREA PLANTING, [0,0.1,...,0.9,1]
for (i in 1:11)
{
  #Loop over release rates [0,4,8,...,36,40]
  for (k in 1:11)
  {
    #Loop over initial resistance levels (in above vector)
    for (ii in 1:9)
    {
      
      
      ###################################################################################################
      # Define Biological Parameters
      ###################################################################################################
      
      #Major Policy Parameters to change - Loop or single run?
      # Loop is i=1:11
      
      .bt.amt = (i-1)/10            # bt.amt = Percentage of crop planted in Bt plants
      #.bt.amt=1
      
      #################################
      # Pre-defining Parameters
      # Parameters to adjust
      
      ###################################################
      #Loop is ii= init.resist.list<-c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.9)
      
      init.resist=init.resist.list[ii]
      #init.resist=0.001        # amount of initial Bt-resistance in wild population (i.e. init.resist)
     
      # release = weekly release ratio of transgenic individuals to wild-type males, r in paper
      # for genotypes 1 - 9
      ##############################
      #Loop or single run?
      .release1 = 4*(k-1)           # Change this value to release tg KKSS
      #.release1 = 36
     
      
      # Theta = Juvenile survival in Bt fields 
      # for genotypes 1 - 9
      .theta1 = 1-(mort.bt*.bt.amt)
      .theta2 = 1-(dom.th*mort.bt*.bt.amt)
      .theta3 = 1
      .theta4 = 1-(mort.bt*.bt.amt)
      .theta5 = 1-(dom.th*mort.bt*.bt.amt)
      .theta6 = 1
      .theta7 = 1-(mort.bt*.bt.amt)
      .theta8 = 1-(dom.th*mort.bt*.bt.amt)
      .theta9 = 1
      
      ## initial genotype frequencies multiplied by the equilibrium value
      # As of 12/2/15: based on expected genotype frequencies in wild population
      # juvenile genotype frequencies
      .j1=equil.j*0
      .j2=equil.j*0
      .j3=equil.j*0
      .j4=equil.j*0
      .j5=equil.j*0
      .j6=equil.j*0
      .j7=equil.j*(1-init.resist)^2
      .j8=equil.j*(2*init.resist*(1-init.resist)) 
      .j9=equil.j*(init.resist^2)
      # female genotype frequencies
      .f1=equil.f*0
      .f2=equil.f*0
      .f3=equil.f*0
      .f4=equil.f*0 
      .f5=equil.f*0
      .f6=equil.f*0
      .f7=equil.f*(1-init.resist)^2
      .f8=equil.f*(2*init.resist*(1-init.resist))
      .f9=equil.f*(init.resist^2)
      # male genotype frequencies
      .m1=equil.m*0 
      .m2=equil.m*0 
      .m3=equil.m*0
      .m4=equil.m*0 
      .m5=equil.m*0
      .m6=equil.m*0
      .m7=equil.m*(1-init.resist)^2
      .m8=equil.m*(2*init.resist*(1-init.resist))
      .m9=equil.m*(init.resist^2)
      
      ###################################################################################################
      # Economics Portion
      ###################################################################################################
      #################################
      
      #Describing parameters
      
      # q =           Refuge [0,1]
      # rho =         Discount Rate [0,1]
      # h =           Dominance factor [0,1]
      # eps =         Generation Scale
      # r =           fitness cost per generation [0,1]
      # mu =          lethality of Bt [0,1]
      # g =           inherent growth rate 
      # K =           Carrying Capacity
      # BtPrice =     Bt seed price, normalized to OUTPUT price
      # InsPrice  =   insecticide price, norm. to OUTPUT price
      # AltSeedPrice= alternate normal seed price, norm to OUTPUT price
      # alpha =       Damage done by juvenile population, norm. by Carrying Capacity
      
      #Describing States
      #phi =          bunch-a-shit
      
      #Y_avg =        average of Bt and Non-Bt area yield, norm. to 1 for "potential yield"
      #Y_bt =         Bt area yield, norm. to 1
      #Y_nonbt =      Non-Bt area yield, norm. to 1
      
      # D =           Juventile population (#)
      # w =           susceptibility [0,1]
      
      # npv =         Net present value of the investment
      
      ##################################
      # Economic Parameters
      #(Most of these aren't used, from Zack's. Keeping for moment)
      .rho=0.07 
#       .h=0 
#       .eps=4 
#       .r=0.2 
#       .mu=0.93 
#       .g=5 
#       .K=1000
#       .w=0.99
#       .D=1000
      
      # Just ball parking.  Maybe not necessary to exactly specify?
      # But higher Bt seed price critical to profit comparison.
#       .BtPrice=0.40 # $/kg
#       .InsPrice= 5 # $/kg
#       .InsAppRate = 2.5 # kg/ha
#       .AltSeedPrice=0.20 # $/kg
#       .OutPrice=0.10 # $/kg - find brocolli prices or price ranges.
#       
      
      
      ################
      #Economics links
      ################
      ### Need to reevaluate which equil used here
      #Total Juveniles in equil.
      .J=equil.j
      #Susceptibility to Bt % in equil.
      .Sus= (.j1 + .j2 + .j4 + .j5 + .j7 + .j8)/.J
      
      #Econ stuff depending on equil value
      .alpha=0.8/15000  #assumption that equil is carrying capacity is WRONG, needs to fix
      #cite: Yi et al. (2015) for only the homozygotes being resistant
      .Y_bt=1-.alpha*(1-.Sus)*.J
      .Y_nonbt=1-(.alpha*.J)
      .Y_avg=(.bt.amt)*.Y_bt+(1-.bt.amt)*.Y_nonbt
      
      .NetBen=(.Y_avg - .Y_nonbt)
      #.NPV=(.NetBen/100)/(1+.rho)^(0)
      
      
      
      ###################################################################################################
      # Define vectors for use in ode function
      ###################################################################################################
      #################################
      # Define length of run - Set at 30 years
      times <- seq(from = 0, to = 50, by = 0.01)
      
      #################################
      # Define the state variables
      state <- c(j1=.j1, j2=.j2, j3=.j3, j4=.j4, j5=.j5, j6=.j6, j7=.j7, j8=.j8, j9=.j9
                 , f1=.f1, f2=.f2, f3=.f3, f4=.f4, f5=.f5, f6=.f6, f7=.f7, f8=.f8, f9=.f9
                 , m1=.m1, m2=.m2, m3=.m3, m4=.m4, m5=.m5, m6=.m6, m7=.m7, m8=.m8, m9=.m9
                 , J=.J, Sus=.Sus
                 , Y_bt=.Y_bt, Y_nonbt=.Y_nonbt, Y_avg=.Y_avg
                 , NetBen=.NetBen
                 #, NPV=.NPV
      )
      
      #################################
      # Define the parameter values
      parameters <- as.list(c(larval.prod=.larval.prod, dens=.dens, strength.dens=.strength.dens
                              , mort.j=.mort.j, mort.f=.mort.f, mort.m=.mort.m, bt.amt=.bt.amt
                              , release1=.release1, release2=.release2, release3=.release3
                              , release4=.release4, release5=.release5, release6=.release6 
                              , release7=.release7, release8=.release8, release9=.release9, M9=.m9
                              , leaky=.leaky, emergence=.emergence
                              , viability1=.viability1, viability2=.viability2, viability3=.viability3
                              , viability4=.viability4, viability5=.viability5, viability6=.viability6 
                              , viability7=.viability7, viability8=.viability8, viability9=.viability9 
                              , fitness1=.fitness1, fitness2=.fitness2, fitness3=.fitness3
                              , fitness4=.fitness4, fitness5=.fitness5, fitness6=.fitness6 
                              , fitness7=.fitness7, fitness8=.fitness8, fitness9=.fitness9 
                              , theta1=.theta1, theta2=.theta2, theta3=.theta3
                              , theta4=.theta4, theta5=.theta5, theta6=.theta6 
                              , theta7=.theta7, theta8=.theta8, theta9=.theta9
                              , rho=.rho
                              #, h=.h, eps=.eps, r=.r, mu=.mu, g=.g, K=.K, BtPrice=.BtPrice
                              #, InsPrice=.InsPrice, AltSeedPrice=.AltSeedPrice, OutPrice=.OutPrice
      ))
      
      ###################################################################################################
      # Define Robert function
      ###################################################################################################
      #################################
      Robert <- function(t, state, parameters){ 
        with(as.list(c(state, parameters)), {
          # where state[1:9] is juvenile density for genotypes 1 - 9
          # where state[10:18] is female density for genotypes 1 - 9
          # where state[19:27] is male density for genotypes 1 - 9
          ### Differential Equations for Genotype 1: KKSS
          dj1 <- 365*((fitness1 * theta1 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[1,1,1]*(state[19]),
                           p.matrix[1,1,2]*(state[20]),
                           p.matrix[1,1,3]*(state[21]),
                           p.matrix[1,1,4]*(state[22]),
                           p.matrix[1,1,5]*(state[23]),
                           p.matrix[1,1,6]*(state[24]),
                           p.matrix[1,1,7]*(state[25]),
                           p.matrix[1,1,8]*(state[26]),
                           p.matrix[1,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[1,2,1]*(state[19]),
                           p.matrix[1,2,2]*(state[20]),
                           p.matrix[1,2,3]*(state[21]),
                           p.matrix[1,2,4]*(state[22]),
                           p.matrix[1,2,5]*(state[23]),
                           p.matrix[1,2,6]*(state[24]),
                           p.matrix[1,2,7]*(state[25]),
                           p.matrix[1,2,8]*(state[26]),
                           p.matrix[1,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[1,3,1]*(state[19]),
                           p.matrix[1,3,2]*(state[20]),
                           p.matrix[1,3,3]*(state[21]),
                           p.matrix[1,3,4]*(state[22]),
                           p.matrix[1,3,5]*(state[23]),
                           p.matrix[1,3,6]*(state[24]),
                           p.matrix[1,3,7]*(state[25]),
                           p.matrix[1,3,8]*(state[26]),
                           p.matrix[1,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[1,4,1]*(state[19]),
                           p.matrix[1,4,2]*(state[20]),
                           p.matrix[1,4,3]*(state[21]),
                           p.matrix[1,4,4]*(state[22]),
                           p.matrix[1,4,5]*(state[23]),
                           p.matrix[1,4,6]*(state[24]),
                           p.matrix[1,4,7]*(state[25]),
                           p.matrix[1,4,8]*(state[26]),
                           p.matrix[1,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[1,5,1]*(state[19]),
                           p.matrix[1,5,2]*(state[20]),
                           p.matrix[1,5,3]*(state[21]),
                           p.matrix[1,5,4]*(state[22]),
                           p.matrix[1,5,5]*(state[23]),
                           p.matrix[1,5,6]*(state[24]),
                           p.matrix[1,5,7]*(state[25]),
                           p.matrix[1,5,8]*(state[26]),
                           p.matrix[1,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[1,6,1]*(state[19]),
                           p.matrix[1,6,2]*(state[20]),
                           p.matrix[1,6,3]*(state[21]),
                           p.matrix[1,6,4]*(state[22]),
                           p.matrix[1,6,5]*(state[23]),
                           p.matrix[1,6,6]*(state[24]),
                           p.matrix[1,6,7]*(state[25]),
                           p.matrix[1,6,8]*(state[26]),
                           p.matrix[1,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[1,7,1]*(state[19]),
                           p.matrix[1,7,2]*(state[20]),
                           p.matrix[1,7,3]*(state[21]),
                           p.matrix[1,7,4]*(state[22]),
                           p.matrix[1,7,5]*(state[23]),
                           p.matrix[1,7,6]*(state[24]),
                           p.matrix[1,7,7]*(state[25]),
                           p.matrix[1,7,8]*(state[26]),
                           p.matrix[1,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[1,8,1]*(state[19]),
                           p.matrix[1,8,2]*(state[20]),
                           p.matrix[1,8,3]*(state[21]),
                           p.matrix[1,8,4]*(state[22]),
                           p.matrix[1,8,5]*(state[23]),
                           p.matrix[1,8,6]*(state[24]),
                           p.matrix[1,8,7]*(state[25]),
                           p.matrix[1,8,8]*(state[26]),
                           p.matrix[1,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[1,9,1]*(state[19]),
                           p.matrix[1,9,2]*(state[20]),
                           p.matrix[1,9,3]*(state[21]),
                           p.matrix[1,9,4]*(state[22]),
                           p.matrix[1,9,5]*(state[23]),
                           p.matrix[1,9,6]*(state[24]),
                           p.matrix[1,9,7]*(state[25]),
                           p.matrix[1,9,8]*(state[26]),
                           p.matrix[1,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[1] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[1]) - (emergence * state[1])
          )
          df1 <- 365*((0.5 * emergence * viability1 * j1) - (mort.f * state[10]) + (release1*leaky*M9/7))
          dm1 <- 365*((0.5 * emergence * j1) - (mort.m * state[19]) + (release1*M9/7))
          
          ### Differential Equations for Genotype 2: KKSs
          dj2 <- 365*((fitness2 * theta2 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[2,1,1]*(state[19]),
                           p.matrix[2,1,2]*(state[20]),
                           p.matrix[2,1,3]*(state[21]),
                           p.matrix[2,1,4]*(state[22]),
                           p.matrix[2,1,5]*(state[23]),
                           p.matrix[2,1,6]*(state[24]),
                           p.matrix[2,1,7]*(state[25]),
                           p.matrix[2,1,8]*(state[26]),
                           p.matrix[2,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[2,2,1]*(state[19]),
                           p.matrix[2,2,2]*(state[20]),
                           p.matrix[2,2,3]*(state[21]),
                           p.matrix[2,2,4]*(state[22]),
                           p.matrix[2,2,5]*(state[23]),
                           p.matrix[2,2,6]*(state[24]),
                           p.matrix[2,2,7]*(state[25]),
                           p.matrix[2,2,8]*(state[26]),
                           p.matrix[2,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[2,3,1]*(state[19]),
                           p.matrix[2,3,2]*(state[20]),
                           p.matrix[2,3,3]*(state[21]),
                           p.matrix[2,3,4]*(state[22]),
                           p.matrix[2,3,5]*(state[23]),
                           p.matrix[2,3,6]*(state[24]),
                           p.matrix[2,3,7]*(state[25]),
                           p.matrix[2,3,8]*(state[26]),
                           p.matrix[2,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[2,4,1]*(state[19]),
                           p.matrix[2,4,2]*(state[20]),
                           p.matrix[2,4,3]*(state[21]),
                           p.matrix[2,4,4]*(state[22]),
                           p.matrix[2,4,5]*(state[23]),
                           p.matrix[2,4,6]*(state[24]),
                           p.matrix[2,4,7]*(state[25]),
                           p.matrix[2,4,8]*(state[26]),
                           p.matrix[2,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[2,5,1]*(state[19]),
                           p.matrix[2,5,2]*(state[20]),
                           p.matrix[2,5,3]*(state[21]),
                           p.matrix[2,5,4]*(state[22]),
                           p.matrix[2,5,5]*(state[23]),
                           p.matrix[2,5,6]*(state[24]),
                           p.matrix[2,5,7]*(state[25]),
                           p.matrix[2,5,8]*(state[26]),
                           p.matrix[2,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[2,6,1]*(state[19]),
                           p.matrix[2,6,2]*(state[20]),
                           p.matrix[2,6,3]*(state[21]),
                           p.matrix[2,6,4]*(state[22]),
                           p.matrix[2,6,5]*(state[23]),
                           p.matrix[2,6,6]*(state[24]),
                           p.matrix[2,6,7]*(state[25]),
                           p.matrix[2,6,8]*(state[26]),
                           p.matrix[2,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[2,7,1]*(state[19]),
                           p.matrix[2,7,2]*(state[20]),
                           p.matrix[2,7,3]*(state[21]),
                           p.matrix[2,7,4]*(state[22]),
                           p.matrix[2,7,5]*(state[23]),
                           p.matrix[2,7,6]*(state[24]),
                           p.matrix[2,7,7]*(state[25]),
                           p.matrix[2,7,8]*(state[26]),
                           p.matrix[2,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[2,8,1]*(state[19]),
                           p.matrix[2,8,2]*(state[20]),
                           p.matrix[2,8,3]*(state[21]),
                           p.matrix[2,8,4]*(state[22]),
                           p.matrix[2,8,5]*(state[23]),
                           p.matrix[2,8,6]*(state[24]),
                           p.matrix[2,8,7]*(state[25]),
                           p.matrix[2,8,8]*(state[26]),
                           p.matrix[2,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[2,9,1]*(state[19]),
                           p.matrix[2,9,2]*(state[20]),
                           p.matrix[2,9,3]*(state[21]),
                           p.matrix[2,9,4]*(state[22]),
                           p.matrix[2,9,5]*(state[23]),
                           p.matrix[2,9,6]*(state[24]),
                           p.matrix[2,9,7]*(state[25]),
                           p.matrix[2,9,8]*(state[26]),
                           p.matrix[2,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[2] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[2]) - (emergence * state[2]))
          df2 <- 365*((0.5 * emergence * viability2 * j2) - (mort.f * state[11]) + (release2*leaky*M9/7))
          dm2 <- 365*((0.5 * emergence * j2) - (mort.m * state[20]) + (release2*M9/7))
          
          ### Differential Equations for Genotype 3: KKss
          dj3 <- 365*((fitness3 * theta3 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[3,1,1]*(state[19]),
                           p.matrix[3,1,2]*(state[20]),
                           p.matrix[3,1,3]*(state[21]),
                           p.matrix[3,1,4]*(state[22]),
                           p.matrix[3,1,5]*(state[23]),
                           p.matrix[3,1,6]*(state[24]),
                           p.matrix[3,1,7]*(state[25]),
                           p.matrix[3,1,8]*(state[26]),
                           p.matrix[3,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[3,2,1]*(state[19]),
                           p.matrix[3,2,2]*(state[20]),
                           p.matrix[3,2,3]*(state[21]),
                           p.matrix[3,2,4]*(state[22]),
                           p.matrix[3,2,5]*(state[23]),
                           p.matrix[3,2,6]*(state[24]),
                           p.matrix[3,2,7]*(state[25]),
                           p.matrix[3,2,8]*(state[26]),
                           p.matrix[3,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[3,3,1]*(state[19]),
                           p.matrix[3,3,2]*(state[20]),
                           p.matrix[3,3,3]*(state[21]),
                           p.matrix[3,3,4]*(state[22]),
                           p.matrix[3,3,5]*(state[23]),
                           p.matrix[3,3,6]*(state[24]),
                           p.matrix[3,3,7]*(state[25]),
                           p.matrix[3,3,8]*(state[26]),
                           p.matrix[3,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[3,4,1]*(state[19]),
                           p.matrix[3,4,2]*(state[20]),
                           p.matrix[3,4,3]*(state[21]),
                           p.matrix[3,4,4]*(state[22]),
                           p.matrix[3,4,5]*(state[23]),
                           p.matrix[3,4,6]*(state[24]),
                           p.matrix[3,4,7]*(state[25]),
                           p.matrix[3,4,8]*(state[26]),
                           p.matrix[3,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[3,5,1]*(state[19]),
                           p.matrix[3,5,2]*(state[20]),
                           p.matrix[3,5,3]*(state[21]),
                           p.matrix[3,5,4]*(state[22]),
                           p.matrix[3,5,5]*(state[23]),
                           p.matrix[3,5,6]*(state[24]),
                           p.matrix[3,5,7]*(state[25]),
                           p.matrix[3,5,8]*(state[26]),
                           p.matrix[3,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[3,6,1]*(state[19]),
                           p.matrix[3,6,2]*(state[20]),
                           p.matrix[3,6,3]*(state[21]),
                           p.matrix[3,6,4]*(state[22]),
                           p.matrix[3,6,5]*(state[23]),
                           p.matrix[3,6,6]*(state[24]),
                           p.matrix[3,6,7]*(state[25]),
                           p.matrix[3,6,8]*(state[26]),
                           p.matrix[3,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[3,7,1]*(state[19]),
                           p.matrix[3,7,2]*(state[20]),
                           p.matrix[3,7,3]*(state[21]),
                           p.matrix[3,7,4]*(state[22]),
                           p.matrix[3,7,5]*(state[23]),
                           p.matrix[3,7,6]*(state[24]),
                           p.matrix[3,7,7]*(state[25]),
                           p.matrix[3,7,8]*(state[26]),
                           p.matrix[3,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[3,8,1]*(state[19]),
                           p.matrix[3,8,2]*(state[20]),
                           p.matrix[3,8,3]*(state[21]),
                           p.matrix[3,8,4]*(state[22]),
                           p.matrix[3,8,5]*(state[23]),
                           p.matrix[3,8,6]*(state[24]),
                           p.matrix[3,8,7]*(state[25]),
                           p.matrix[3,8,8]*(state[26]),
                           p.matrix[3,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[3,9,1]*(state[19]),
                           p.matrix[3,9,2]*(state[20]),
                           p.matrix[3,9,3]*(state[21]),
                           p.matrix[3,9,4]*(state[22]),
                           p.matrix[3,9,5]*(state[23]),
                           p.matrix[3,9,6]*(state[24]),
                           p.matrix[3,9,7]*(state[25]),
                           p.matrix[3,9,8]*(state[26]),
                           p.matrix[3,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[3] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[3]) - (emergence * state[3]))
          df3 <- 365*((0.5 * emergence * viability3 * j3) - (mort.f * state[12]) + (release3*leaky*M9/7))
          dm3 <- 365*((0.5 * emergence * j3) - (mort.m * state[21]) + (release3*M9/7))
          
          ### Differential Equations for Genotype 4: KkSS
          dj4 <- 365*((fitness4 * theta4 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[4,1,1]*(state[19]),
                           p.matrix[4,1,2]*(state[20]),
                           p.matrix[4,1,3]*(state[21]),
                           p.matrix[4,1,4]*(state[22]),
                           p.matrix[4,1,5]*(state[23]),
                           p.matrix[4,1,6]*(state[24]),
                           p.matrix[4,1,7]*(state[25]),
                           p.matrix[4,1,8]*(state[26]),
                           p.matrix[4,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[4,2,1]*(state[19]),
                           p.matrix[4,2,2]*(state[20]),
                           p.matrix[4,2,3]*(state[21]),
                           p.matrix[4,2,4]*(state[22]),
                           p.matrix[4,2,5]*(state[23]),
                           p.matrix[4,2,6]*(state[24]),
                           p.matrix[4,2,7]*(state[25]),
                           p.matrix[4,2,8]*(state[26]),
                           p.matrix[4,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[4,3,1]*(state[19]),
                           p.matrix[4,3,2]*(state[20]),
                           p.matrix[4,3,3]*(state[21]),
                           p.matrix[4,3,4]*(state[22]),
                           p.matrix[4,3,5]*(state[23]),
                           p.matrix[4,3,6]*(state[24]),
                           p.matrix[4,3,7]*(state[25]),
                           p.matrix[4,3,8]*(state[26]),
                           p.matrix[4,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[4,4,1]*(state[19]),
                           p.matrix[4,4,2]*(state[20]),
                           p.matrix[4,4,3]*(state[21]),
                           p.matrix[4,4,4]*(state[22]),
                           p.matrix[4,4,5]*(state[23]),
                           p.matrix[4,4,6]*(state[24]),
                           p.matrix[4,4,7]*(state[25]),
                           p.matrix[4,4,8]*(state[26]),
                           p.matrix[4,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[4,5,1]*(state[19]),
                           p.matrix[4,5,2]*(state[20]),
                           p.matrix[4,5,3]*(state[21]),
                           p.matrix[4,5,4]*(state[22]),
                           p.matrix[4,5,5]*(state[23]),
                           p.matrix[4,5,6]*(state[24]),
                           p.matrix[4,5,7]*(state[25]),
                           p.matrix[4,5,8]*(state[26]),
                           p.matrix[4,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[4,6,1]*(state[19]),
                           p.matrix[4,6,2]*(state[20]),
                           p.matrix[4,6,3]*(state[21]),
                           p.matrix[4,6,4]*(state[22]),
                           p.matrix[4,6,5]*(state[23]),
                           p.matrix[4,6,6]*(state[24]),
                           p.matrix[4,6,7]*(state[25]),
                           p.matrix[4,6,8]*(state[26]),
                           p.matrix[4,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[4,7,1]*(state[19]),
                           p.matrix[4,7,2]*(state[20]),
                           p.matrix[4,7,3]*(state[21]),
                           p.matrix[4,7,4]*(state[22]),
                           p.matrix[4,7,5]*(state[23]),
                           p.matrix[4,7,6]*(state[24]),
                           p.matrix[4,7,7]*(state[25]),
                           p.matrix[4,7,8]*(state[26]),
                           p.matrix[4,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[4,8,1]*(state[19]),
                           p.matrix[4,8,2]*(state[20]),
                           p.matrix[4,8,3]*(state[21]),
                           p.matrix[4,8,4]*(state[22]),
                           p.matrix[4,8,5]*(state[23]),
                           p.matrix[4,8,6]*(state[24]),
                           p.matrix[4,8,7]*(state[25]),
                           p.matrix[4,8,8]*(state[26]),
                           p.matrix[4,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[4,9,1]*(state[19]),
                           p.matrix[4,9,2]*(state[20]),
                           p.matrix[4,9,3]*(state[21]),
                           p.matrix[4,9,4]*(state[22]),
                           p.matrix[4,9,5]*(state[23]),
                           p.matrix[4,9,6]*(state[24]),
                           p.matrix[4,9,7]*(state[25]),
                           p.matrix[4,9,8]*(state[26]),
                           p.matrix[4,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[4] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[4]) - (emergence * state[4]))
          df4 <- 365*((0.5 * emergence * viability4 * j4) - (mort.f * state[13]) + (release4*leaky*M9/7))
          dm4 <- 365*((0.5 * emergence * j4) - (mort.m * state[22]) + (release4*M9/7))
          
          ### Differential Equations for Genotype 5: KkSs
          dj5 <- 365*((fitness5 * theta5 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[5,1,1]*(state[19]),
                           p.matrix[5,1,2]*(state[20]),
                           p.matrix[5,1,3]*(state[21]),
                           p.matrix[5,1,4]*(state[22]),
                           p.matrix[5,1,5]*(state[23]),
                           p.matrix[5,1,6]*(state[24]),
                           p.matrix[5,1,7]*(state[25]),
                           p.matrix[5,1,8]*(state[26]),
                           p.matrix[5,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[5,2,1]*(state[19]),
                           p.matrix[5,2,2]*(state[20]),
                           p.matrix[5,2,3]*(state[21]),
                           p.matrix[5,2,4]*(state[22]),
                           p.matrix[5,2,5]*(state[23]),
                           p.matrix[5,2,6]*(state[24]),
                           p.matrix[5,2,7]*(state[25]),
                           p.matrix[5,2,8]*(state[26]),
                           p.matrix[5,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[5,3,1]*(state[19]),
                           p.matrix[5,3,2]*(state[20]),
                           p.matrix[5,3,3]*(state[21]),
                           p.matrix[5,3,4]*(state[22]),
                           p.matrix[5,3,5]*(state[23]),
                           p.matrix[5,3,6]*(state[24]),
                           p.matrix[5,3,7]*(state[25]),
                           p.matrix[5,3,8]*(state[26]),
                           p.matrix[5,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[5,4,1]*(state[19]),
                           p.matrix[5,4,2]*(state[20]),
                           p.matrix[5,4,3]*(state[21]),
                           p.matrix[5,4,4]*(state[22]),
                           p.matrix[5,4,5]*(state[23]),
                           p.matrix[5,4,6]*(state[24]),
                           p.matrix[5,4,7]*(state[25]),
                           p.matrix[5,4,8]*(state[26]),
                           p.matrix[5,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[5,5,1]*(state[19]),
                           p.matrix[5,5,2]*(state[20]),
                           p.matrix[5,5,3]*(state[21]),
                           p.matrix[5,5,4]*(state[22]),
                           p.matrix[5,5,5]*(state[23]),
                           p.matrix[5,5,6]*(state[24]),
                           p.matrix[5,5,7]*(state[25]),
                           p.matrix[5,5,8]*(state[26]),
                           p.matrix[5,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[5,6,1]*(state[19]),
                           p.matrix[5,6,2]*(state[20]),
                           p.matrix[5,6,3]*(state[21]),
                           p.matrix[5,6,4]*(state[22]),
                           p.matrix[5,6,5]*(state[23]),
                           p.matrix[5,6,6]*(state[24]),
                           p.matrix[5,6,7]*(state[25]),
                           p.matrix[5,6,8]*(state[26]),
                           p.matrix[5,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[5,7,1]*(state[19]),
                           p.matrix[5,7,2]*(state[20]),
                           p.matrix[5,7,3]*(state[21]),
                           p.matrix[5,7,4]*(state[22]),
                           p.matrix[5,7,5]*(state[23]),
                           p.matrix[5,7,6]*(state[24]),
                           p.matrix[5,7,7]*(state[25]),
                           p.matrix[5,7,8]*(state[26]),
                           p.matrix[5,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[5,8,1]*(state[19]),
                           p.matrix[5,8,2]*(state[20]),
                           p.matrix[5,8,3]*(state[21]),
                           p.matrix[5,8,4]*(state[22]),
                           p.matrix[5,8,5]*(state[23]),
                           p.matrix[5,8,6]*(state[24]),
                           p.matrix[5,8,7]*(state[25]),
                           p.matrix[5,8,8]*(state[26]),
                           p.matrix[5,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[5,9,1]*(state[19]),
                           p.matrix[5,9,2]*(state[20]),
                           p.matrix[5,9,3]*(state[21]),
                           p.matrix[5,9,4]*(state[22]),
                           p.matrix[5,9,5]*(state[23]),
                           p.matrix[5,9,6]*(state[24]),
                           p.matrix[5,9,7]*(state[25]),
                           p.matrix[5,9,8]*(state[26]),
                           p.matrix[5,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[5] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[5]) - (emergence * state[5]))
          df5 <- 365*((0.5 * emergence * viability5 * j5) - (mort.f * state[14]) + (release5*leaky*M9/7))
          dm5 <- 365*((0.5 * emergence * j5) - (mort.m * state[23]) + (release5*M9/7))
          
          ### Differential Equations for Genotype 6: Kkss
          dj6 <- 365*((fitness6 * theta6 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[6,1,1]*(state[19]),
                           p.matrix[6,1,2]*(state[20]),
                           p.matrix[6,1,3]*(state[21]),
                           p.matrix[6,1,4]*(state[22]),
                           p.matrix[6,1,5]*(state[23]),
                           p.matrix[6,1,6]*(state[24]),
                           p.matrix[6,1,7]*(state[25]),
                           p.matrix[6,1,8]*(state[26]),
                           p.matrix[6,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[6,2,1]*(state[19]),
                           p.matrix[6,2,2]*(state[20]),
                           p.matrix[6,2,3]*(state[21]),
                           p.matrix[6,2,4]*(state[22]),
                           p.matrix[6,2,5]*(state[23]),
                           p.matrix[6,2,6]*(state[24]),
                           p.matrix[6,2,7]*(state[25]),
                           p.matrix[6,2,8]*(state[26]),
                           p.matrix[6,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[6,3,1]*(state[19]),
                           p.matrix[6,3,2]*(state[20]),
                           p.matrix[6,3,3]*(state[21]),
                           p.matrix[6,3,4]*(state[22]),
                           p.matrix[6,3,5]*(state[23]),
                           p.matrix[6,3,6]*(state[24]),
                           p.matrix[6,3,7]*(state[25]),
                           p.matrix[6,3,8]*(state[26]),
                           p.matrix[6,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[6,4,1]*(state[19]),
                           p.matrix[6,4,2]*(state[20]),
                           p.matrix[6,4,3]*(state[21]),
                           p.matrix[6,4,4]*(state[22]),
                           p.matrix[6,4,5]*(state[23]),
                           p.matrix[6,4,6]*(state[24]),
                           p.matrix[6,4,7]*(state[25]),
                           p.matrix[6,4,8]*(state[26]),
                           p.matrix[6,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[6,5,1]*(state[19]),
                           p.matrix[6,5,2]*(state[20]),
                           p.matrix[6,5,3]*(state[21]),
                           p.matrix[6,5,4]*(state[22]),
                           p.matrix[6,5,5]*(state[23]),
                           p.matrix[6,5,6]*(state[24]),
                           p.matrix[6,5,7]*(state[25]),
                           p.matrix[6,5,8]*(state[26]),
                           p.matrix[6,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[6,6,1]*(state[19]),
                           p.matrix[6,6,2]*(state[20]),
                           p.matrix[6,6,3]*(state[21]),
                           p.matrix[6,6,4]*(state[22]),
                           p.matrix[6,6,5]*(state[23]),
                           p.matrix[6,6,6]*(state[24]),
                           p.matrix[6,6,7]*(state[25]),
                           p.matrix[6,6,8]*(state[26]),
                           p.matrix[6,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[6,7,1]*(state[19]),
                           p.matrix[6,7,2]*(state[20]),
                           p.matrix[6,7,3]*(state[21]),
                           p.matrix[6,7,4]*(state[22]),
                           p.matrix[6,7,5]*(state[23]),
                           p.matrix[6,7,6]*(state[24]),
                           p.matrix[6,7,7]*(state[25]),
                           p.matrix[6,7,8]*(state[26]),
                           p.matrix[6,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[6,8,1]*(state[19]),
                           p.matrix[6,8,2]*(state[20]),
                           p.matrix[6,8,3]*(state[21]),
                           p.matrix[6,8,4]*(state[22]),
                           p.matrix[6,8,5]*(state[23]),
                           p.matrix[6,8,6]*(state[24]),
                           p.matrix[6,8,7]*(state[25]),
                           p.matrix[6,8,8]*(state[26]),
                           p.matrix[6,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[6,9,1]*(state[19]),
                           p.matrix[6,9,2]*(state[20]),
                           p.matrix[6,9,3]*(state[21]),
                           p.matrix[6,9,4]*(state[22]),
                           p.matrix[6,9,5]*(state[23]),
                           p.matrix[6,9,6]*(state[24]),
                           p.matrix[6,9,7]*(state[25]),
                           p.matrix[6,9,8]*(state[26]),
                           p.matrix[6,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[6] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[6]) - (emergence * state[6]))
          df6 <- 365*((0.5 * emergence * viability6 * j6) - (mort.f * state[15]) + (release6*leaky*M9/7))
          dm6 <- 365*((0.5 * emergence * j6) - (mort.m * state[24]) + (release6*M9/7))
          
          ### Differential Equations for Genotype 7: kkSS
          dj7 <- 365*((fitness7 * theta7 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[7,1,1]*(state[19]),
                           p.matrix[7,1,2]*(state[20]),
                           p.matrix[7,1,3]*(state[21]),
                           p.matrix[7,1,4]*(state[22]),
                           p.matrix[7,1,5]*(state[23]),
                           p.matrix[7,1,6]*(state[24]),
                           p.matrix[7,1,7]*(state[25]),
                           p.matrix[7,1,8]*(state[26]),
                           p.matrix[7,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[7,2,1]*(state[19]),
                           p.matrix[7,2,2]*(state[20]),
                           p.matrix[7,2,3]*(state[21]),
                           p.matrix[7,2,4]*(state[22]),
                           p.matrix[7,2,5]*(state[23]),
                           p.matrix[7,2,6]*(state[24]),
                           p.matrix[7,2,7]*(state[25]),
                           p.matrix[7,2,8]*(state[26]),
                           p.matrix[7,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[7,3,1]*(state[19]),
                           p.matrix[7,3,2]*(state[20]),
                           p.matrix[7,3,3]*(state[21]),
                           p.matrix[7,3,4]*(state[22]),
                           p.matrix[7,3,5]*(state[23]),
                           p.matrix[7,3,6]*(state[24]),
                           p.matrix[7,3,7]*(state[25]),
                           p.matrix[7,3,8]*(state[26]),
                           p.matrix[7,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[7,4,1]*(state[19]),
                           p.matrix[7,4,2]*(state[20]),
                           p.matrix[7,4,3]*(state[21]),
                           p.matrix[7,4,4]*(state[22]),
                           p.matrix[7,4,5]*(state[23]),
                           p.matrix[7,4,6]*(state[24]),
                           p.matrix[7,4,7]*(state[25]),
                           p.matrix[7,4,8]*(state[26]),
                           p.matrix[7,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[7,5,1]*(state[19]),
                           p.matrix[7,5,2]*(state[20]),
                           p.matrix[7,5,3]*(state[21]),
                           p.matrix[7,5,4]*(state[22]),
                           p.matrix[7,5,5]*(state[23]),
                           p.matrix[7,5,6]*(state[24]),
                           p.matrix[7,5,7]*(state[25]),
                           p.matrix[7,5,8]*(state[26]),
                           p.matrix[7,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[7,6,1]*(state[19]),
                           p.matrix[7,6,2]*(state[20]),
                           p.matrix[7,6,3]*(state[21]),
                           p.matrix[7,6,4]*(state[22]),
                           p.matrix[7,6,5]*(state[23]),
                           p.matrix[7,6,6]*(state[24]),
                           p.matrix[7,6,7]*(state[25]),
                           p.matrix[7,6,8]*(state[26]),
                           p.matrix[7,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[7,7,1]*(state[19]),
                           p.matrix[7,7,2]*(state[20]),
                           p.matrix[7,7,3]*(state[21]),
                           p.matrix[7,7,4]*(state[22]),
                           p.matrix[7,7,5]*(state[23]),
                           p.matrix[7,7,6]*(state[24]),
                           p.matrix[7,7,7]*(state[25]),
                           p.matrix[7,7,8]*(state[26]),
                           p.matrix[7,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[7,8,1]*(state[19]),
                           p.matrix[7,8,2]*(state[20]),
                           p.matrix[7,8,3]*(state[21]),
                           p.matrix[7,8,4]*(state[22]),
                           p.matrix[7,8,5]*(state[23]),
                           p.matrix[7,8,6]*(state[24]),
                           p.matrix[7,8,7]*(state[25]),
                           p.matrix[7,8,8]*(state[26]),
                           p.matrix[7,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[7,9,1]*(state[19]),
                           p.matrix[7,9,2]*(state[20]),
                           p.matrix[7,9,3]*(state[21]),
                           p.matrix[7,9,4]*(state[22]),
                           p.matrix[7,9,5]*(state[23]),
                           p.matrix[7,9,6]*(state[24]),
                           p.matrix[7,9,7]*(state[25]),
                           p.matrix[7,9,8]*(state[26]),
                           p.matrix[7,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[7] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[7]) - (emergence * state[7]))
          df7 <- 365*((0.5 * emergence * viability7 * j7) - (mort.f * state[16]) + (release7*leaky*M9/7))
          dm7 <- 365*((0.5 * emergence * j7) - (mort.m * state[25]) + (release7*M9/7))
          
          ### Differential Equations for Genotype 8: kkSs
          dj8 <- 365*((fitness8 * theta8 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[8,1,1]*(state[19]),
                           p.matrix[8,1,2]*(state[20]),
                           p.matrix[8,1,3]*(state[21]),
                           p.matrix[8,1,4]*(state[22]),
                           p.matrix[8,1,5]*(state[23]),
                           p.matrix[8,1,6]*(state[24]),
                           p.matrix[8,1,7]*(state[25]),
                           p.matrix[8,1,8]*(state[26]),
                           p.matrix[8,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[8,2,1]*(state[19]),
                           p.matrix[8,2,2]*(state[20]),
                           p.matrix[8,2,3]*(state[21]),
                           p.matrix[8,2,4]*(state[22]),
                           p.matrix[8,2,5]*(state[23]),
                           p.matrix[8,2,6]*(state[24]),
                           p.matrix[8,2,7]*(state[25]),
                           p.matrix[8,2,8]*(state[26]),
                           p.matrix[8,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[8,3,1]*(state[19]),
                           p.matrix[8,3,2]*(state[20]),
                           p.matrix[8,3,3]*(state[21]),
                           p.matrix[8,3,4]*(state[22]),
                           p.matrix[8,3,5]*(state[23]),
                           p.matrix[8,3,6]*(state[24]),
                           p.matrix[8,3,7]*(state[25]),
                           p.matrix[8,3,8]*(state[26]),
                           p.matrix[8,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[8,4,1]*(state[19]),
                           p.matrix[8,4,2]*(state[20]),
                           p.matrix[8,4,3]*(state[21]),
                           p.matrix[8,4,4]*(state[22]),
                           p.matrix[8,4,5]*(state[23]),
                           p.matrix[8,4,6]*(state[24]),
                           p.matrix[8,4,7]*(state[25]),
                           p.matrix[8,4,8]*(state[26]),
                           p.matrix[8,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[8,5,1]*(state[19]),
                           p.matrix[8,5,2]*(state[20]),
                           p.matrix[8,5,3]*(state[21]),
                           p.matrix[8,5,4]*(state[22]),
                           p.matrix[8,5,5]*(state[23]),
                           p.matrix[8,5,6]*(state[24]),
                           p.matrix[8,5,7]*(state[25]),
                           p.matrix[8,5,8]*(state[26]),
                           p.matrix[8,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[8,6,1]*(state[19]),
                           p.matrix[8,6,2]*(state[20]),
                           p.matrix[8,6,3]*(state[21]),
                           p.matrix[8,6,4]*(state[22]),
                           p.matrix[8,6,5]*(state[23]),
                           p.matrix[8,6,6]*(state[24]),
                           p.matrix[8,6,7]*(state[25]),
                           p.matrix[8,6,8]*(state[26]),
                           p.matrix[8,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[8,7,1]*(state[19]),
                           p.matrix[8,7,2]*(state[20]),
                           p.matrix[8,7,3]*(state[21]),
                           p.matrix[8,7,4]*(state[22]),
                           p.matrix[8,7,5]*(state[23]),
                           p.matrix[8,7,6]*(state[24]),
                           p.matrix[8,7,7]*(state[25]),
                           p.matrix[8,7,8]*(state[26]),
                           p.matrix[8,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[8,8,1]*(state[19]),
                           p.matrix[8,8,2]*(state[20]),
                           p.matrix[8,8,3]*(state[21]),
                           p.matrix[8,8,4]*(state[22]),
                           p.matrix[8,8,5]*(state[23]),
                           p.matrix[8,8,6]*(state[24]),
                           p.matrix[8,8,7]*(state[25]),
                           p.matrix[8,8,8]*(state[26]),
                           p.matrix[8,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[8,9,1]*(state[19]),
                           p.matrix[8,9,2]*(state[20]),
                           p.matrix[8,9,3]*(state[21]),
                           p.matrix[8,9,4]*(state[22]),
                           p.matrix[8,9,5]*(state[23]),
                           p.matrix[8,9,6]*(state[24]),
                           p.matrix[8,9,7]*(state[25]),
                           p.matrix[8,9,8]*(state[26]),
                           p.matrix[8,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[8] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[8]) - (emergence * state[8]))
          df8 <- 365*((0.5 * emergence * viability8 * j8) - (mort.f * state[17]) + (release8*leaky*M9/7))
          dm8 <- 365*((0.5 * emergence * j8) - (mort.m * state[26]) + (release8*M9/7))
          
          ### Differential Equations for Genotype 9: kkss (aka: wild-type)
          dj9 <- 365*((fitness9 * theta9 * larval.prod *
                     sum(
                       state[10]* 
                         sum(
                           p.matrix[9,1,1]*(state[19]),
                           p.matrix[9,1,2]*(state[20]),
                           p.matrix[9,1,3]*(state[21]),
                           p.matrix[9,1,4]*(state[22]),
                           p.matrix[9,1,5]*(state[23]),
                           p.matrix[9,1,6]*(state[24]),
                           p.matrix[9,1,7]*(state[25]),
                           p.matrix[9,1,8]*(state[26]),
                           p.matrix[9,1,9]*(state[27])
                         ),
                       
                       state[11]*
                         sum(
                           p.matrix[9,2,1]*(state[19]),
                           p.matrix[9,2,2]*(state[20]),
                           p.matrix[9,2,3]*(state[21]),
                           p.matrix[9,2,4]*(state[22]),
                           p.matrix[9,2,5]*(state[23]),
                           p.matrix[9,2,6]*(state[24]),
                           p.matrix[9,2,7]*(state[25]),
                           p.matrix[9,2,8]*(state[26]),
                           p.matrix[9,2,9]*(state[27])
                         ),
                       
                       state[12]*
                         sum(
                           p.matrix[9,3,1]*(state[19]),
                           p.matrix[9,3,2]*(state[20]),
                           p.matrix[9,3,3]*(state[21]),
                           p.matrix[9,3,4]*(state[22]),
                           p.matrix[9,3,5]*(state[23]),
                           p.matrix[9,3,6]*(state[24]),
                           p.matrix[9,3,7]*(state[25]),
                           p.matrix[9,3,8]*(state[26]),
                           p.matrix[9,3,9]*(state[27])
                         ),
                       
                       state[13]*
                         sum(
                           p.matrix[9,4,1]*(state[19]),
                           p.matrix[9,4,2]*(state[20]),
                           p.matrix[9,4,3]*(state[21]),
                           p.matrix[9,4,4]*(state[22]),
                           p.matrix[9,4,5]*(state[23]),
                           p.matrix[9,4,6]*(state[24]),
                           p.matrix[9,4,7]*(state[25]),
                           p.matrix[9,4,8]*(state[26]),
                           p.matrix[9,4,9]*(state[27])
                         ),
                       
                       state[14]*
                         sum(
                           p.matrix[9,5,1]*(state[19]),
                           p.matrix[9,5,2]*(state[20]),
                           p.matrix[9,5,3]*(state[21]),
                           p.matrix[9,5,4]*(state[22]),
                           p.matrix[9,5,5]*(state[23]),
                           p.matrix[9,5,6]*(state[24]),
                           p.matrix[9,5,7]*(state[25]),
                           p.matrix[9,5,8]*(state[26]),
                           p.matrix[9,5,9]*(state[27])
                         ),
                       
                       state[15]*
                         sum(
                           p.matrix[9,6,1]*(state[19]),
                           p.matrix[9,6,2]*(state[20]),
                           p.matrix[9,6,3]*(state[21]),
                           p.matrix[9,6,4]*(state[22]),
                           p.matrix[9,6,5]*(state[23]),
                           p.matrix[9,6,6]*(state[24]),
                           p.matrix[9,6,7]*(state[25]),
                           p.matrix[9,6,8]*(state[26]),
                           p.matrix[9,6,9]*(state[27])
                         ),
                       
                       state[16]*
                         sum(
                           p.matrix[9,7,1]*(state[19]),
                           p.matrix[9,7,2]*(state[20]),
                           p.matrix[9,7,3]*(state[21]),
                           p.matrix[9,7,4]*(state[22]),
                           p.matrix[9,7,5]*(state[23]),
                           p.matrix[9,7,6]*(state[24]),
                           p.matrix[9,7,7]*(state[25]),
                           p.matrix[9,7,8]*(state[26]),
                           p.matrix[9,7,9]*(state[27])
                         ),
                       
                       state[17]*
                         sum(
                           p.matrix[9,8,1]*(state[19]),
                           p.matrix[9,8,2]*(state[20]),
                           p.matrix[9,8,3]*(state[21]),
                           p.matrix[9,8,4]*(state[22]),
                           p.matrix[9,8,5]*(state[23]),
                           p.matrix[9,8,6]*(state[24]),
                           p.matrix[9,8,7]*(state[25]),
                           p.matrix[9,8,8]*(state[26]),
                           p.matrix[9,8,9]*(state[27])
                         ),
                       
                       state[18]*
                         sum(
                           p.matrix[9,9,1]*(state[19]),
                           p.matrix[9,9,2]*(state[20]),
                           p.matrix[9,9,3]*(state[21]),
                           p.matrix[9,9,4]*(state[22]),
                           p.matrix[9,9,5]*(state[23]),
                           p.matrix[9,9,6]*(state[24]),
                           p.matrix[9,9,7]*(state[25]),
                           p.matrix[9,9,8]*(state[26]),
                           p.matrix[9,9,9]*(state[27])
                         )
                     )/sum(state[19:27])
          )
          - state[9] * (dens * sum(state[1:9]))^(strength.dens - 1)
          - (mort.j * state[9]) - (emergence * state[9]))
          df9 <- 365*((0.5 * emergence * viability9 * j9) - (mort.f * state[18]) + (release9*leaky*M9/7))
          dm9 <- 365*((0.5 * emergence * j9) - (mort.m * state[27]) + (release9*M9/7))
          
          ##Economics links
          #Damage per insect (was originally "alpha=0.8/(carrying capacity)" in Zack's)
          alpha=.alpha
          
          # What would be the dominance factor (h in Zack's) equivalent?
          
          Y_bt=1-alpha*(1-Sus)*J
          
          #Need to incorporate a demand function that only plants Bt when 
          #more profitable than non-Bt 
          #i.e. Demand= 1-q if more profitable, and less if not. 
          # Likely need cont. for compl vs sub question.
          #Also need to decompose alpha into more interpretable per Juv. damage. 
          #Use mkt rejection paper and cross-ref rejection rates with field pop.s
          
          Y_nonbt=1-(alpha*J)
          Y_avg=(bt.amt)*Y_bt+(1-bt.amt)*Y_nonbt
          
          NetBen=(Y_avg - Y_nonbt)
          #Zack said just to do Y_avg
          
          #This may be too specific with prices for below.
          # ProfitGain=(Y_avg*OutPrice-q*BtPrice-(1-q)*AltSeedPrice)-(Y_nonbt*OutPrice-AltSeedPrice)
          #^^ clearly wrong, but fix tomorrow.  Focus on Net Yield Ben. for now.
          
          
          
          ####Originals for combining with genetics section    
          #DE for Total juveniles in pop causing crop damage
          dJ <- dj1 + dj2 + dj3 + dj4 + dj5 + dj6 + dj7 + dj8 + dj9
          # DE for % with susceptibility to Bt
          dSus <- (dj1 + dj2 + dj4 + dj5 + dj7 + dj8)/J - ((j1+j2+j4+j5+j7+j8)/(J^2))*dJ
          
          #Total differentiation by hand, I think on target but potential errors. Check with Wolfram Alpha
          dY_bt <- (-1)*alpha*(1-Sus)*dJ + alpha*J*dSus
          dY_nonbt <- (-1)*alpha*dJ
          dY_avg <- (bt.amt)*dY_bt+(1-bt.amt)*dY_nonbt
          dNetBen <- dY_avg - Y_nonbt
          # dNPV <- ((dY_avg - dY_nonbt)/100)/((1+rho)^t)
          
          # Return results
          list(c(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9, df1, df2, df3, df4, 
                 df5, df6, df7, df8, df9, dm1, dm2, dm3, dm4, dm5, dm6, dm7, dm8, dm9
                 ,dJ, dSus, dY_bt, dY_nonbt, dY_avg
                    , dNetBen
                 #   , dNPV
          ))
        }) # End with statement
      } # end function
      
      #################################
      # Solve for differential equations in Robert function
      out <- ode(y = state, times = times, func = Robert, parms = parameters)
      head(out)
      
      # #####
      #Cumulative NPV for Bt
      npv.cum=rep(0,length(out[,1]))
      
      for (j in 1:length(out[,1]))
      {
        npv.cum[j]=(out[j,33]/100)/((1+.rho)^(out[j,1]))
        
      }
      sum(npv.cum)
      
      # Plot NPV vs Time
   #   par(oma = c(0, 0, 3, 0))
    #  plot(cumsum(npv.cum),xlab="time",ylab="Cumulative NPV")
      
      ##########################################################
      ##Loop series or single run?
      npv.cum.list[i,k,ii]=sum(npv.cum)
      
    } #End loop for initial resistance level
  } #End loop for release number
} #End loop for refuge area (q)


npv.cum.list

######################################################
# # Plotting Genotype Results Individually
#par(mar=c(4,3,3,3), mgp=c(2,1,0))
# plot(out, xlab = "Time")

# ######################################################
# # Plotting Full Results
# par(mar=c(4,3,3,3), mgp=c(2,1,0))
# 
# t_out=out[,1]
# 
# KKSS_num=out[,2]
# KKSs_num=out[,3]
# KKss_num=out[,4]
# 
# KkSS_num=out[,5]
# KkSs_num=out[,6]
# Kkss_num=out[,7]
# 
# kkSS_num=out[,8]
# kkSs_num=out[,9]
# kkss_num=out[,10]
# 
# N_total=rowSums(out[,8:10])
# 
# # Plot Number of individuals from each genotype over time
# pdf("Release40_kkSS_Bt.amt1.pdf")
# plot(x=NULL, y=NULL, xlim=c(0,700), ylim=c(0,10000), 
#      xlab="Time", ylab="Population Size", type='n', 
#      main="Releasing Wild-Type Bt Susceptibles", 
#      sub="Release Ratio = 40, Cost.s = 0.03, Bt.amt = 1, init.resist = 0.06")
# 
# lines(t_out, N_total, type='l', col="black")
# lines(t_out, kkSS_num, type='l', col="red")
# lines(t_out, kkSs_num, type='l', col="green")
# lines(t_out, kkss_num, type='l', col="blue")
# legend("right", c("Total","kkSS","kkSs","kkss"), lty=c(1,1), 
#        col=c("black","red","green","blue"), title="Genotypes")
# dev.off()

# ######################################################
# # Create contour plot of npv.cum.list
# par(mar = c(2, 4, 2, 2))
# image2D(npv.cum.list, rasterImage = TRUE, contour = list(lwd = 2, col = jet.col(11)))
# title(main = "NPV Contour Plots", sub = NULL, xlab = "Refuge", ylab = "Test")

