import numpy as np
import matplotlib.pyplot as plt

# define the function
def f(x):
    return (np.power(x[0],4)-2*np.power(x[0],2)+x[0])+(np.power(x[1],4)-2*np.power(x[1],2)+x[1])

# define constants to be used
n = 2
iterations = 30
m = 100
radius = 2.0
sigma = 5.0
alpha = 1.0

# define initial conditions
x = np.zeros((iterations,2))
fx = np.zeros(iterations)
x[0] = [1.1,1.2]
fx[0] = f(x[0])

# draw contour plot
xlist = np.linspace(-2.0, 2.0, 100) # Create 1-D arrays for x,y dimensions
ylist = np.linspace(-2.0, 2.0, 100)
X,Y = np.meshgrid(xlist, ylist) # Create 2-D grid xlist,ylist values
Z = f(np.array([X,Y]))
plt.contour(X, Y, Z, [-5.0,-4.0,-3.9,-3.8,-3.7,-3.6,-3.5,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0], colors = 'b', linestyles = 'solid')

#iteration
for k in range(0,iterations-1):
    u = np.zeros((m,n))
    count = 0
    while count < m:
        unew = np.zeros(n)
        for i in range(0,n):
            unew[i] = np.random.normal(x[k,i],sigma)
        distance = np.sum(np.abs(unew-x[k])**2,axis=-1)**(1./2)
        if distance < radius:
            u[count] = unew
            count += 1

    fu = np.array([f(u_) for u_ in u])
    j = np.argmin(fu)

    # accept the new variable if function value gets better
    if fu[j] < fx[k]:
        x[k+1] = u[j]
        fx[k+1] = f(x[k+1])
    else:
        p = np.zeros(m)
        for i in range(0,m):
            p[i] = np.exp(alpha*(fx[k]-fu[i]))
        S = np.sum(p)
        p = p/S
        xi = np.random.rand()
        for i in range(0,m):
            if xi < np.sum(p[:i]):
                x[k+1] = u[i]
                fx[k+1] = f(x[k+1])

bestX = x[np.argmin(fx)]
bestF = np.min(fx)

# print the function values changing
print fx
print "Best variable result:"
print bestX
print "Best function value:"
print bestF

annotation = "f(%3.8f,%3.8f)=%3.8f" % (bestX[0],bestX[1],bestF)
# plot
plt.plot(x[:,0],x[:,1])
plt.title("Simulated Annealing Simple 2D Example \n with f = x^4-2x^2+x+y^4-2y^2+y")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(bestX[0],bestX[1],'^')
plt.annotate(annotation, (bestX[0],bestX[1]))
plt.show()
