######## exercise 9 #####
# Dan and Drew




#1

import pandas
from plotnine import *
import scipy
pondat = pandas.read_csv("ponzr1.csv")

#barplot  for mean counts per mutation
grouped= pondat.groupby(["mutation"]).mean().reset_index() #mean counts by mutation
print grouped

grouped.columns = ['mutation', 'mean_counts']
p= (ggplot(data=grouped)
    + aes(x='mutation', y= 'mean_counts',fill= 'mutation')
    + geom_bar(stat = "identity")
    + theme_classic()
    )
print p


subset_WT = pondat.query('mutation=="WT"')
subset_M124K = pondat.query('mutation=="M124K"')
subset_V456D = pondat.query('mutation=="V456D"')
subset_I213N = pondat.query('mutation=="I213N"')






#2





#3
