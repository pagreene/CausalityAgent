import re
import os
import logging
import sqlite3
import numpy as np


class CausalityAgent:

    def __init__(self):
        # Load CA_target database
        correlation_db_file = './resources/PNNL-ovarian-correlations.db'
        if os.path.isfile(correlation_db_file):
            self.correlation_db = sqlite3.connect(correlation_db_file)

        else:
            self.correlation_db = sqlite3.connect(correlation_db_file)
            self.populate_correlation_db()
            # self.correlation_db = None;

        causality_db_file = './resources/causality.db'
        if os.path.isfile(causality_db_file):
            self.causality_db = sqlite3.connect(causality_db_file)

        else:
            open(causality_db_file, 'w')
            self.causality_db = sqlite3.connect(causality_db_file)
            self.populate_causality_db()
            # self.causality_db = None;

        self.geneInd = 0

    def __del__(self):
        self.correlation_db.close()
        self.causality_db.close()

    def populate_causality_db(self):
        causality_file = open("./resources/causative-data-centric.sif", 'r')

        with self.causality_db:
            cur = self.causality_db.cursor()
            try:
                cur.execute(
                    "CREATE TABLE Causality(Id1 TEXT, PSite1 TEXT, Id2 TEXT, PSite2 TEXT, Rel TEXT, UriStr TEXT)")
            except:
                pass

            for line in causality_file:

                vals = line.split('\t')
                id_str1 = vals[0].upper().split('-')
                id1 = id_str1[0]
                if len(id_str1) > 1:
                    p_site1 = id_str1[1]
                else:
                    p_site1 = ' '

                id_str2 = vals[2].upper().split('-')
                id2 = id_str2[0]

                if len(id_str2) > 1:
                    p_site2 = id_str2[1]
                else:
                    p_site2 = ' '

                rel = vals[1]

                uri_arr = []
                if vals[3]:
                    uri_arr = vals[3].split(" ")

                if len(uri_arr) == 0:
                    uri_arr = [vals[3]]

                uri_str = ""
                for uri in uri_arr:
                    uri_str = uri_str + "uri= " + uri + "&"


                print((id1, p_site1, id2, p_site2, rel, uri_str))
                cur.execute("INSERT INTO Causality VALUES(?, ?, ?, ?, ?, ?)", (id1, p_site1, id2, p_site2, rel, uri_str))

        causality_file.close()


    def populate_correlation_db(self):
        pnnl_file = open("./resources/PNNL-ovarian-correlations.txt", 'r')

        with self.correlation_db:
            cur = self.correlation_db.cursor()
            try:
                cur.execute(
                    "CREATE TABLE Correlations(Id1 TEXT, PSite1 TEXT, Id2 TEXT, PSite2 TEXT, Corr REAL, PVal REAL)")
            except:
                pass

            for line in pnnl_file:
                if line.find('/') > -1:  # incorrectly formatted strings
                    continue
                vals = line.split('\t')
                id_str1 = vals[0].upper().split('-')
                id1 = id_str1[0]
                if len(id_str1) > 1:
                    p_site1 = id_str1[1]
                else:
                    p_site1 = ' '

                    id_str2 = vals[1].upper().split('-')
                id2 = id_str2[0]

                if len(id_str2) > 1:
                    p_site2 = id_str2[1]
                else:
                    p_site2 = ' '

                corr = float(vals[2].rstrip('\n'))

                p_val = float(vals[3].rstrip('\n'))

                cur.execute("INSERT INTO Correlations VALUES(?, ?, ?, ?, ?, ?)", (id1, p_site1, id2, p_site2, corr, p_val))

        pnnl_file.close()

    # Find the n'th highest correlation
    def find_next_correlation(self, gene):

        with self.correlation_db:
            cur = self.correlation_db.cursor()
            rows = cur.execute("SELECT * FROM Correlations WHERE Id1 = ?  ORDER BY ABS(Corr) DESC LIMIT ?",
                               (gene, self.geneInd + 1)).fetchall()
            # for row in rows:
            #     print row[0], row[1], row[2], row[3], row[4]

            row_cnt = len(rows)
            if row_cnt > 0:
                self.geneInd = self.geneInd + 1
                return rows[row_cnt-1]

    def find_causality(self, gene1, gene2):
        with self.causality_db:
            cur = self.causality_db.cursor()
            rows = cur.execute("SELECT * FROM Causality WHERE Id1 = ? AND  Id2 = ?  OR  Id1 = ? AND Id2 = ?",
                               (gene1, gene2, gene2, gene1)).fetchall()

            return rows



agent = CausalityAgent()
# agent.populate_correlation_db()
# agent.populate_causality_db()
corr = agent.find_next_correlation('AKT1')
print(corr[2])
corr = agent.find_next_correlation('AKT1')
print(corr[2])
causality = agent.find_causality('MAPK1', 'JUND')
print(causality[0][4])
