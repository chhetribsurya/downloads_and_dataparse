
"""Index of duplicates items in a python list"""

source = "ABABDBAAEDSBQEWBAFLSAFB"

from collections import defaultdict

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() 
                            if len(locs)>1)

for dup in sorted(list_duplicates(source)):
    print dup

Prints:
('A', [0, 2, 6, 7, 16, 20])
('B', [1, 3, 5, 11, 15, 22])
('D', [4, 9])
('E', [8, 13])
('F', [17, 21])
('S', [10, 19])

#If you want to do repeated testing for various keys against the same source, you can use functools.partial to create a new function variable, using a "partially complete" argument list, that is, specifying the seq, but omitting the item to search for:

from functools import partial
dups_in_source = partial(list_duplicates_of, source)

for c in "ABDEFS":
    print c, dups_in_source(c)

Prints:
A [0, 2, 6, 7, 16, 20]
B [1, 3, 5, 11, 15, 22]
D [4, 9]
E [8, 13]
F [17, 21]
S [10, 19]


dups = collections.defaultdict(list)
for i, e in enumerate(L):
  dups[e].append(i)
for k, v in sorted(dups.iteritems()):
  if len(v) >= 2:
    print '%s: %r' % (k, v)


# def duplicates(lst, item):
# 	return [i for i, x in enumerate(lst) if x == item]

# duplicates(List, "A")
# [0, 2]



""" Tested and proven, use this concept with defaultdict """

In [1]: source = list("abcabababcdddeeefff")

In [3]: from collections import defaultdict

In [4]: dup_dict = defaultdict(list)

In [5]: dup_dict
Out[5]: defaultdict(list, {})

In [7]: for i,item in enumerate(source):
            dup_dict[item].append(i)

In [8]: dup_dict
Out[8]: 
defaultdict(list,
            {'a': [0, 3, 5, 7],
             'b': [1, 4, 6, 8],
             'c': [2, 9],
             'd': [10, 11, 12],
             'e': [13, 14, 15],
             'f': [16, 17, 18]})


In [12]: dup_dict["a"]
Out[12]: [0, 3, 5, 7]

In [13]: dup_dict.values()
Out[13]: 
[[0, 3, 5, 7],
 [],
 [2, 9],
 [1, 4, 6, 8],
 [13, 14, 15],
 [10, 11, 12],
 [],
 [16, 17, 18]]



In [14]: dups_dict = {}

In [15]: for i,item in enumerate(source):
   ....:     if item in dups_dict:
   ....:         dups_dict[item].append(i)
   ....:     else:
   ....:         dups_dict[item] = []
   ....:         

In [16]: dups_dict
Out[16]: 
{'a': [3, 5, 7],
 'b': [4, 6, 8],
 'c': [9],
 'd': [11, 12],
 'e': [14, 15],
 'f': [17, 18]}

In [17]: dups_dict = {}

In [18]: for i,item in enumerate(source):
    if item in dups_dict:
        dups_dict[item].append(i)
    else:
        dups_dict[item] = [i]
   ....:         

In [19]: dups_dict 
Out[19]: 
{'a': [0, 3, 5, 7],
 'b': [1, 4, 6, 8],
 'c': [2, 9],
 'd': [10, 11, 12],
 'e': [13, 14, 15],
 'f': [16, 17, 18]}





#2 – Apply Function

#It is one of the commonly used functions for playing with data and creating new variables. Apply returns some value after passing each row/column of a data frame with some function. The function can be both default or user-defined. For instance, here it can be used to find the #missing values in each row and column.

#Create a new function:
def num_missing(x):
  return sum(x.isnull())

#Applying per column:
print "Missing values per column:"
print data.apply(num_missing, axis=0) #axis=0 defines that function is to be applied on each column

#Applying per row:
print "\nMissing values per row:"
print data.apply(num_missing, axis=1).head() #axis=1 defines that function is to be applied on each row




#10 – Cut function for binning

"""Sometimes numerical values make more sense if clustered together. For example, if we’re trying to model traffic (#cars on road) with time of the day (minutes). The exact minute of an hour might not be that relevant for predicting traffic as compared to actual period of the day like “Morning”, “Afternoon”, “Evening”, “Night”, “Late Night”. Modeling traffic this way will be more intuitive and will avoid overfitting.

Here we define a simple function which can be re-used for binning any variable fairly easily."""

#Binning:
def binning(col, cut_points, labels=None):
  #Define min and max values:
  minval = col.min()
  maxval = col.max()

  #create list by adding min and max to cut_points
  break_points = [minval] + cut_points + [maxval]

  #if no labels provided, use default labels 0 ... (n-1)
  if not labels:
    labels = range(len(cut_points)+1)

  #Binning using cut function of pandas
  colBin = pd.cut(col,bins=break_points,labels=labels,include_lowest=True)
  return colBin

#Binning age:
cut_points = [90,140,190]
labels = ["low","medium","high","very high"]
data["LoanAmount_Bin"] = binning(data["LoanAmount"], cut_points, labels)
print pd.value_counts(data["LoanAmount_Bin"], sort=False)




#3 – Imputing missing files

"""‘fillna()’ does it in one go. It is used for updating missing values with the overall mean/mode/median of the column. Let’s impute the ‘Gender’, ‘Married’ and ‘Self_Employed’ columns with their respective modes."""

#First we import a function to determine the mode
from scipy.stats import mode
mode(data['Gender'])
Output: ModeResult(mode=array([‘Male’], dtype=object), count=array([489]))

"""This returns both mode and count. Remember that mode can be an array as there can be multiple values with high frequency. We will take the first one by default always using:"""

mode(data['Gender']).mode[0]
4.2 - mode2

"""Now we can fill the missing values and check using technique #2."""

#Impute the values:
data['Gender'].fillna(mode(data['Gender']).mode[0], inplace=True)
data['Married'].fillna(mode(data['Married']).mode[0], inplace=True)
data['Self_Employed'].fillna(mode(data['Self_Employed']).mode[0], inplace=True)

#Now check the #missing values again to confirm:
print data.apply(num_missing, axis=0)
