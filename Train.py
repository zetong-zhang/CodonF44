# Author: Zhang Zetong
# Date: 2024-08-17
import pickle
from Bio import SeqIO
from random import shuffle
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


def coding(s):
    p = []
    s, le = str(s).upper(), len(s)
    c = {b: 0 for b in "AGCT"}

    for b in s:
        c[b] += 1

    x = (c['A'] + c['G'] - c['C'] - c['T']) / le
    y = (c['A'] + c['C'] - c['G'] - c['T']) / le
    z = (c['A'] + c['T'] - c['G'] - c['C']) / le

    p += [x, y, z]

    c = {i: {b: 0 for b in "AGCT"} for i in range(3)}

    for i in range(le):
        c[i % 3][s[i]] += 1

    for i in range(3):
        su = sum(c[i].values())
        x = (c[i]['A'] + c[i]['G'] - c[i]['C'] - c[i]['T']) / su
        y = (c[i]['A'] + c[i]['C'] - c[i]['G'] - c[i]['T']) / su
        z = (c[i]['A'] + c[i]['T'] - c[i]['G'] - c[i]['C']) / su
        p += [x, y, z]

    c = {a: {b: 0 for b in "AGCT"} for a in "AGCT"}

    for i in range(le - 1):
        c[s[i]][s[i + 1]] += 1

    for a in "AGCT":
        su = sum(c[a].values())
        if su == 0:
            p += [0, 0, 0]
        else:
            x = (c[a]['A'] + c[a]['G'] - c[a]['C'] - c[a]['T']) / su
            y = (c[a]['A'] + c[a]['C'] - c[a]['G'] - c[a]['T']) / su
            z = (c[a]['A'] + c[a]['T'] - c[a]['G'] - c[a]['C']) / su
            p += [x, y, z]

    return p


records = SeqIO.parse(handle='F44.CDS.fa', format='fasta')
X, Y = [], []

for record in records:
    X.append(coding(record.seq))
    Y.append(1)
    lst = list(record.seq)
    shuffle(lst)
    X.append(coding(''.join(lst)))
    Y.append(0)

model = SVC()
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)
print(round(roc_auc_score(y_test, y_pred), 4))

model.fit(X, Y)

with open("CodonF44_Model.pkl", 'wb') as file:
    pickle.dump(model, file)
