import csv
import os

from copy import deepcopy

def write_final(outfile, named):
    '''
        Joins graph meta-data contained in graphdata.csv on
        graph names resulting from tests written to outfile.
    '''

    tmp_file = 'tmp.csv'

    if named:
        # read graph meta-data from the provided csv
        graph_names = dict()
        with open("graphdata.csv", "r") as graphdatacsv:
            graphreader = csv.reader(graphdatacsv,delimiter=',')

            for row in graphreader:
                if row[0] not in graph_names:
                    graph_names[row[0].replace(" ","")] = deepcopy(row)

        # inner join the meta-data with the results on graph name
        # and write those contents to the outfile provided
        written = False
        with open(outfile, "r") as results,\
            open(tmp_file, "w") as final:

            resultsreader = csv.reader(results,delimiter=',')
            finalwriter = csv.writer(final,delimiter=',')

            for row in resultsreader:
                name = row[0].split('.')[0].replace(" ", "")
                if name in graph_names:
                    all_in_one = graph_names[name] + row
                    finalwriter.writerow(all_in_one)
                else:
                    print(f'graph {name} missing')

    else:
        with open(outfile, "r") as results,\
            open(tmp_file, "w") as final:

            resultsreader = csv.reader(results,delimiter=',')
            finalwriter = csv.writer(final,delimiter=',')

            for row in resultsreader:
                finalwriter.writerow(row)

    os.rename(tmp_file, outfile)