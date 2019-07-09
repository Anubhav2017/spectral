import matplotlib.pyplot as plt
x=[]*none
y=[]*none
i=0
f=open("badFaceParameters","r+")
for line in f:
    line=f.readline()
    for var in line.split():
        print(var)
        if(i%2  == 0):
            x.append(var)
        else:
            y.append(var)

plt.scatter(x,y)
plt.show()