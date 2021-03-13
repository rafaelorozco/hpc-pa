from numpy import random

n = 100

file1 = open("sample-constants"+str(n)+".txt","w") 

file1.write(str(n))
file1.write("\n") 
for i in range(n):
	file1.write(str(random.randint(100)) )
	file1.write("\n") 
file1.close()
