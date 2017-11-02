######## exercise 9 #####
# Dan and Drew
import scipy 
import numpy
import os

# 1
# import libs
import pandas
from plotnine import *
import numpy
from scipy.optimize import minimize
from scipy.stats import norm
import scipy.stats

# load dataframe
pondat = pandas.read_csv("ponzr1.csv")

# barplot  for mean counts per mutation
grouped= pondat.groupby(["mutation"]).mean().reset_index() # mean counts by mutation
grouped.columns = ['mutation', 'mean_counts']
p= (ggplot(data=grouped)
    + aes(x='mutation', y= 'mean_counts',fill= 'mutation')
    + geom_bar(stat = "identity")
    + theme_classic()
    )
print p

# format and subset dataframes
pondat.columns = ['mutation', 'y'] # rename col
pondat['x'] = [0 if ele  == "WT" else 1 for ele in pondat["mutation"]] # assign row values depending on mutation

pondat.set_index("mutation", inplace=True) # reset index

M124K = pondat.loc[['M124K','WT']] # subset for each mutation
V456D = pondat.loc[['V456D','WT']]
I213N = pondat.loc[['I213N','WT']]

# null function
def null(p, obs):
    B0 = p[0]
    sigma = p[1]
    expected = B0
    nll = -1*norm(expected, sigma).logpdf(obs.y).sum()
    return nll

# alt function
def alter (p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    expected =B0+B1*obs.x
    nll = -1*norm(expected, sigma).logpdf(obs.y).sum()
    return nll

# initial guess
initialGuess = numpy.array([1,1,1])

# mk124k
fitNull = minimize(null, initialGuess, method="Nelder-Mead", options={'disp':True}, args=M124K)
fitAlter = minimize(alter, initialGuess, method="Nelder-Mead", options={'disp':True}, args=M124K)
D = 2*(fitNull.fun-fitAlter.fun)
print "MK124K"
print 1 - scipy.stats.chi2.cdf(x=D, df=1) #insig

# v456D
fitNull = minimize(null, initialGuess, method="Nelder-Mead", options={'disp':True}, args=V456D)
fitAlter = minimize(alter, initialGuess, method="Nelder-Mead", options={'disp':True}, args=V456D)
D = 2*(fitNull.fun-fitAlter.fun)
print "V456D"
print 1 - scipy.stats.chi2.cdf(x=D, df=1) #sig treatment

# I212N
fitNull = minimize(null, initialGuess, method="Nelder-Mead", options={'disp':True}, args=I213N)
fitAlter = minimize(alter, initialGuess, method="Nelder-Mead", options={'disp':True}, args=I213N)
D = 2*(fitNull.fun-fitAlter.fun)
print "I212N"
print 1 - scipy.stats.chi2.cdf(x=D, df=1) # insig






#2
bakturriuh = pandas.read_csv("MmarinumGrowth.csv")

def nullbac(p, obs):
    B0 = p[0]
    sigma = p[1]
    expected = mumax*(S/S+Ks)
    nll = -1*norm(expected, sigma).logpdf(obs.y).sum()
    return nll



#3

leaf = pandas.read_csv('leafDecomp.csv') # load data
d = (ggplot(data=leaf)
     + aes(x= "Ms", y= "decomp")
     + geom_point()
     )
print d

# simple functions
def simple (p, obs):
    B0 = p[0]
    sigma = p[1]
    expected = B0
    nll = -1*norm(expected, sigma).logpdf(obs.decomp).sum()
    return nll

# alt1 function
def complex (p, obs):
    B0 = p[0]
    B1 = p[1]
    sigma = p[2]
    expected =B0+B1*obs.Ms
    nll = -1*norm(expected, sigma).logpdf(obs.decomp).sum()
    return nll

# more complex alt
def morecomplex (p, obs):
    B0 = p[0]
    B1 = p[1]
    B2 = p[2]
    sigma = p[3]
    expected =B0+(B1*(obs.Ms))+(B2*(obs.Ms*obs.Ms))
    nll = -1*norm(expected, sigma).logpdf(obs.decomp).sum()
    return nll

initialGuess = numpy.array([1,1,1,1])
comGuess = numpy.array([200,10,-.2,1])

# null to linear
fitNull = minimize(simple, initialGuess, method="Nelder-Mead", options={'disp':True}, args=leaf)
fitAlter = minimize(complex, initialGuess, method="Nelder-Mead", options={'disp':True}, args=leaf)
print "null or simple model values"
print fitNull #these values are right compared to Stuart's answers
print "complex or linear model values"
print fitAlter # these values are right compared to Stuart's answers
D = 2*(fitNull.fun-fitAlter.fun)
print "simple to linear model = sig!"
print 1 - scipy.stats.chi2.cdf(x=D, df=1)

# linear to quadratic
fitNull = minimize(complex, initialGuess, method="Nelder-Mead", options={'disp':True}, args=leaf)
fitAlter = minimize(morecomplex, comGuess, method="Nelder-Mead", options={'disp':True}, args=leaf)
print "complex/linear model values"
print fitNull
print "quadratic model values"
print fitAlter
print "linear model to quadratic = sig!"
print 1 - scipy.stats.chi2.cdf(x=D, df=1)


# null to quadratic
fitNull = minimize(simple, initialGuess, method="Nelder-Mead", options={'disp':True}, args=leaf)
fitAlter = minimize(morecomplex, comGuess, method="Nelder-Mead", options={'disp':True}, args=leaf)
print "null values"
print fitNull
print "quadratic model values"
print fitAlter
print "simple model to quadratic = sig!"
print 1 - scipy.stats.chi2.cdf(x=D, df=2)


## according to our calcs, we should use the hump shaped model!!