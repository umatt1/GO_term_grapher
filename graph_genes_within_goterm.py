import argparse  # Specify arguments of the GO Term to find genes for and the file to find logfoldchanges for
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('GO_Term', help='What GO Term to find genes within')
parser.add_argument('File1', help='What file located in the same directory to find logfold changes of genes in')
parser.add_argument('File2', help='What file located in the same directory to find logfold changes of genes in')
args = parser.parse_args()
File1 = args.File1
File2 = args.File2
GO_TERM = args.GO_Term


# this is the content of GOTERM_INVESTIGATOR in the form of a function. Finds genes in go term in data
def GOTERM_INVESTIGATOR(file, cool_terms):
    df = pd.read_csv(file, sep='\t')
    df = df.drop(labels=['baseMean', 'Gene', 'pvalue', 'SE log2FC', 'Wald Stats'], axis=1)
    # drop terms in dataframe that aren't involved in GO term
    df = df[df['GeneID'].isin(cool_terms)]

    print(df)
    return df


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
    m_value, b_value = np.polyfit(df['logfoldchanges_set1'], df['logfoldchanges_set2'], 1)

    return m_value, b_value


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
        if (df['log2FoldChange_x'][digits] > 0 and df['log2FoldChange_y'][digits] < 0) and (
                df['padj_x'][digits] < .05 or df['padj_y'][digits] < .05):
            color = 'red'
        elif (df['log2FoldChange_x'][digits] < 0) and (df['log2FoldChange_y'][digits] > 0) and (
                df['padj_x'][digits] < .05 or df['padj_y'][digits] < .05):
            color = 'blue'
        else:
            color = "grey"
        plt.scatter(x=df['log2FoldChange_x'][digits], y=df['log2FoldChange_y'][digits], color=color)  # color by p-adj
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
    plt.text(x=minimum-.05, y=minimum-.05, s="y = x", fontsize=5)

    point1y = minimum * m + b
    point2y = maximum * m + b
    plt.plot([minimum, maximum], [point1y, point2y])
    plt.text(x=minimum-.05,y=point1y-.05,s=f"y = {str(m)[:5]}x + {str(b)[:7]}", fontsize=5)
    plt.plot([-(abs(minimum)), +abs(maximum)], [-abs(minimum), +abs(maximum)])

    plt.xlabel(file2)
    plt.ylabel(file1)
    plt.xlim(minimum-.1, maximum+.1)
    plt.ylim(minimum-.1, maximum+.1)
    plt.grid()

    definitions = open("drosophila_names.txt")
    for row in definitions: # label what the gene actually does
        lineArr = row.rstrip().split("\t")
        if lineArr[0] == term:
            term += (" (" + lineArr[1] + ") ")
            break
        else:
            continue

    plt.title(f"Expression differences between {file1} and {file2} for {term}")
    plt.show()


cool_terms = set()  # Specify a list to add the genes we care about

genelist = open('gene_associations_foripage.txt')  # open list of GO:Terms and genes to compare

for line in genelist:  # look for genes in the GO term
    linearr = line.rstrip().split("\t")
    if GO_TERM in linearr:  # if the gene is involved in the GO term we care about, let's remember it
        cool_terms.add(linearr[0])
    else:  # otherwise, let's not
        continue
genelist.close()

print("\nFinding genes, logfoldchanges, padj values of " + GO_TERM + " in " + File1 + "... \n")
df1 = GOTERM_INVESTIGATOR(File1, cool_terms)

print("\nFinding genes, logfoldchanges, padj values of " + GO_TERM + " in " + File2 + "... \n")
df2 = GOTERM_INVESTIGATOR(File2, cool_terms)

m, b = linereg(File1, File2)
print("\nLine of regression: y = " + str(m) + "x + " + str(b) + "\n")

scatplotter(df1, df2, m, b, File1, File2, GO_TERM)
