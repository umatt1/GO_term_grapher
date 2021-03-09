import argparse  # Specify arguments of the GO Term to find genes for and the file to find logfoldchanges for
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('GO_Term', help='What GO terms should be looked for?')
parser.add_argument('File1', help='What file located in the same directory to find logfold changes of genes in')
parser.add_argument('File2', help='What file located in the same directory to find logfold changes of genes in')
args = parser.parse_args()
File1 = args.File1
File2 = args.File2
GO_TERM = args.GO_Term
input = open(GO_TERM)


# takes averages of data for dataframe1 and adds it to dataframe2
def TERM_AVERAGER(dataframe1, dataframe2, goTerm):
    logchanges = []
    padjusted = []

    for i in dataframe1["log2FoldChange"]:  # add to lists
        logchanges.append(i)
    for i in dataframe1["padj"]:
        padjusted.append(i)

    # calculate averages from lists
    averagelogchange = 0
    averagepadjusted = 0
    for i in logchanges:
        averagelogchange += i
    for i in padjusted:
        averagepadjusted += i
    averagelogchange = averagelogchange / len(logchanges)
    averagepadjusted = averagepadjusted / len(padjusted)

    dataframe3 = pd.DataFrame.from_dict(
        {"GeneID": goTerm, "log2FoldChange": averagelogchange, "padj": averagepadjusted})
    dataframe2.append(dataframe3)


# this is the content of linereg in the form of a function. Generates line of regression of overall data
def linereg(file1, file2):
    input1 = open(file1)  # Specify the input file
    input2 = open(file2)

    logfoldchanges1 = []
    logfoldchanges2 = []
    genes1 = []
    genes2 = []
    digits = 0

    for line in input1:  # look line by line through the sorted list of genes within GO term of interest
        lineArr = line.rstrip().split("\t")
        if digits == 0:
            digits = 1
            continue
        else:
            genes1.append(lineArr[0])
            try:
                logfoldchanges1.append(float(lineArr[2]))
            except ValueError:
                logfoldchanges1.append(pd.NaT)
    digits = 0
    for line in input2:  # look line by line through the sorted list of genes within GO term of interest
        lineArr = line.rstrip().split("\t")
        if digits == 0:
            digits = 1
            continue
        else:
            genes2.append(lineArr[0])
            try:
                logfoldchanges2.append(float(lineArr[2]))
            except ValueError:
                logfoldchanges2.append(pd.NaT)

    # convert to pandas dataframe
    set1 = pd.DataFrame.from_dict({'gene': genes1, 'logfoldchanges_set1': logfoldchanges1})
    set1 = set1.dropna()
    set2 = pd.DataFrame.from_dict({'gene': genes2, 'logfoldchanges_set2': logfoldchanges2})
    set2 = set2.dropna()

    # inner merge dataframes
    df = set1.merge(set2, on='gene')
    df = df.astype({'logfoldchanges_set1': 'float', 'logfoldchanges_set2': 'float'})
    m, b = np.polyfit(df['logfoldchanges_set1'], df['logfoldchanges_set2'], 1)

    return m, b


# this is the content of scatplotter in the form of a function. Generates scatter plot
def scatplotter(dataframe1, dataframe2, m, b, file1, file2, term):
    print(dataframe1)
    print(dataframe2)

    print('\nMerging Dataframes...')

    dataframe3 = dataframe1.merge(dataframe2, on='GeneID')

    print('\nDropping NAs...')

    df = dataframe3.dropna()
    df.reset_index(inplace=True)

    print(df)

    # plot points
    df = df.astype({'log2FoldChange_x': 'float', 'log2FoldChange_y': 'float', 'padj_y': 'float', 'padj_x': 'float'})
    digits = 0
    for gene in df['GeneID']:
        if df['padj_x'][digits] > .2 and df['padj_y'][digits] > .2:
            color = 'red'
        elif df['padj_x'][digits] > .2 or df['padj_y'][digits] > .2:
            color = 'blue'
        else:
            color = "grey"
        plt.scatter(x=df['log2FoldChange_x'][digits], y=df['log2FoldChange_y'][digits], color=color)
        i = random.uniform(-.01, .01)
        if color == "grey":
            color = "black"
        plt.text(x=df['log2FoldChange_x'][digits], y=df['log2FoldChange_y'][digits]+i, s=gene, fontsize=5, color=color)
        digits += 1

    # draw standard
    maximum = 0.0
    minimum = 0.0
    for number in df['log2FoldChange_x']:
        if number > maximum:
            maximum = number
        if number < minimum:
            minimum = number
    for number in df['log2FoldChange_y']:
        if number > maximum:
            maximum = number
        if number < minimum:
            minimum = number

    point1y = minimum * m + b
    point2y = maximum * m + b
    plt.plot([minimum, maximum], [point1y, point2y])
    plt.plot([-(abs(minimum)), +abs(maximum)], [-abs(minimum), +abs(maximum)])
    plt.text(x=minimum-.05,y=point1y-.05,s=f"y = {str(m)[:5]}x + {str(b)[:7]}", fontsize=5)
    plt.text(x=minimum-.05, y=minimum-.05, s="y = x", fontsize=5)


    plt.xlabel(file2)
    plt.ylabel(file1)
    plt.xlim(minimum-.1,maximum+.1)
    plt.ylim(minimum-.1,maximum+.1)
    plt.grid()
    plt.title("Expression differences between " + file1 + " and " + file2 + " for " + term)
    plt.show()

# open the list of specified GO terms, add them to a list
termlist = set()
digits = 0
for line in input:  # look for genes in the GO term
    if digits == 0:
        digits = 1
        continue
    else:
        termlist.add(line[:-1])

# create dataframes to add to
columns = {"GeneID": [], "log2FoldChange": [], "padj": []}
df1averages = pd.DataFrame.from_dict(columns)
df2averages = pd.DataFrame.from_dict(columns)

dictionary = {}  # dictionary to add the terms of GO term to
genelist = open('gene_associations_foripage.txt')  # open list of GO:Terms and genes to compare
for line in genelist:
    linearr = line.rstrip().split("\t")
    for term in termlist:
        for item in linearr:
            if term == item:
                dictionary.update({linearr[0]: term})

logchanges1 = {}
padj1 = {}
for term in termlist:
    logchanges1.update({term: []})
    padj1.update({term: []})

for line in open(File1):
    linearr = line.rstrip().split("\t")
    if linearr[0] in dictionary.keys():
        try:
            logchanges1[dictionary[linearr[0]]].append(float(linearr[2]))
        except ValueError:
            logchanges1[dictionary[linearr[0]]].append((linearr[2]))

        try:
            padj1[dictionary[linearr[0]]].append(float(linearr[5]))
        except ValueError:
            logchanges1[dictionary[linearr[0]]].append((linearr[2]))

for term in logchanges1:
    total = 0.0
    items = 0
    for number in logchanges1[term]:
        try:
            total += float(number)
        except ValueError:
            total += 0.0
        items += 1
    try: logchanges1.update({term: (total / items)})
    except ZeroDivisionError: logchanges1.update({term: 0})

for term in padj1:
    total = 0
    significantNums = 0
    for number in padj1[term]:
        total += 1
        if float(number) < .05:
            significantNums += 1
    try: padj1.update({term: (significantNums / total)})
    except ZeroDivisionError: padj1.update({term: 0})


# ---- REPEAT FOR NEXT FILE

dictionary = {}  # dictionary to add the terms of GO term to
genelist = open('gene_associations_foripage.txt')  # open list of GO:Terms and genes to compare
for line in genelist:
    linearr = line.rstrip().split("\t")
    for term in termlist:
        for item in linearr:
            if term == item:
                dictionary.update({linearr[0]: term})

logchanges2 = {}
padj2 = {}
for term in termlist:
    logchanges2.update({term: []})
    padj2.update({term: []})

for line in open(File2):
    linearr = line.rstrip().split("\t")
    if linearr[0] in dictionary.keys():
        try:
            logchanges2[dictionary[linearr[0]]].append(float(linearr[2]))
        except ValueError:
            logchanges2[dictionary[linearr[0]]].append((linearr[2]))

        try:
            padj2[dictionary[linearr[0]]].append(float(linearr[5]))
        except ValueError:
            logchanges2[dictionary[linearr[0]]].append((linearr[2]))

for term in logchanges2:
    total = 0.0
    items = 0
    for number in logchanges2[term]:
        try:
            total += float(number)
        except ValueError:
            total += 0.0
        items += 1
    try: logchanges2.update({term: (total / items)})
    except ZeroDivisionError: logchanges2.update({term: 0})

for term in padj2:
    total = 0
    significantNums = 0
    for number in padj2[term]:
        total += 1
        if float(number) < .05:
            significantNums += 1
    try: padj2.update({term: (significantNums / total)})
    except ZeroDivisionError: padj2.update({term: 0})

columns = ["GeneID", "log2FoldChange", "padj"]


# add data to df1averages and df2averages dataframes
for term in logchanges1:
    dataframe = pd.DataFrame.from_dict({0: [term, logchanges1[term], padj1[term]]}, orient='index', columns=columns)
    df1averages = df1averages.append(dataframe)
for term in logchanges2:
    dataframe = pd.DataFrame.from_dict({0: [term, logchanges2[term], padj2[term]]}, orient='index', columns=columns)
    df2averages = df2averages.append(dataframe)

print(df1averages)
print(df2averages)

m, b = linereg(File1, File2)
print("\nLine of regression: y = " + str(m) + "x + " + str(b) + "\n")

scatplotter(df1averages, df2averages, m, b, File1, File2, "selected GO terms")

