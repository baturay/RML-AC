# transform DATASET3.csv to the format we need

f = open("DATASET3.csv")
lines = f.readlines()
g = open("D3", "w")
for i, line in enumerate(lines):
    toks = line.split(",")

    cl = ""
    if i < 400:
        cl = "1"
    elif i < 600:
        cl = "2"
    elif i < 800:
        cl = "3"
    elif i < 1000:
        cl = "4"
    elif i < 1200:
        cl = "5"
    elif i < 1400:
        cl = "6"
    elif i < 1600:
        cl = "4"
    elif i < 2200:
        cl = "6"

    g.write(cl + "," + ",".join(toks[:2]) + "\n")

g.close()
