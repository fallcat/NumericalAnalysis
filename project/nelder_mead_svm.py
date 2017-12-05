import numpy as np
import matplotlib.pyplot as plt

# define the function
def f(x):
    return np.power(x[0],2) + np.power(x[1],2)

# define constants to be used
alpha = 1.0
gamma = 2.0
beta = 0.5
epsilon = 0.01

# initialize x array and other variables
x = np.array([[100.0,150.0],[120.0,-90.0],[-110.0,-130.0]])
fx = np.array([f(x_) for x_ in x])
fxsort = fx.argsort()
fx = fx[fxsort]
x = x[fxsort]
n = 2
count = 0

# draw contour plot
xlist = np.linspace(-1.0, 1.0, 100) # Create 1-D arrays for x,y dimensions
ylist = np.linspace(-1.0, 1.0, 100)
X,Y = np.meshgrid(xlist, ylist) # Create 2-D grid xlist,ylist values
Z = f(np.array([X,Y]))
plt.contour(X, Y, Z, [10000.0, 20000.0, 30000.0, 40000.0], colors = 'b', linestyles = 'solid')

# iteration
# reflection
while True:
    count += 1
    xnew = []
    xbar = 1/float(n)*np.sum(x[:n-1,:],axis=0) # centroid
    xr = (1+alpha)*xbar - alpha*x[n]
    fxr = f(xr)
    if fx[0] <= fxr <= fx[n-1]:
        xnew = xr #reflection_accepted = true
        print "reflection"
    elif fxr <= fx[0]:
        # expansion
        xe = gamma*xr + (1-gamma)*xbar
        if f(xe) < fx[0]:
            xnew = xe # expansion_accepted = true
            print "expansion"
        else:
            xnew = xr # reflection_accepted = true
            print "reflection"
    else: # fx[n-1] <= fxr:
        # contraction
        if fx[n] <= fxr:
            # internal contraction
            xc = beta*fx[n]+(1-beta)*xbar
        else:
            # external contraction
            xc = beta*xr + (1-beta)*xbar
        if f(xc) < fx[n-1]:
            xnew = xc
            print "contraction"
        else: # both reflection vertex and contraction vertex are rejected
            # shrinkage
            for i in range(1,n):
                x[i] = (x[i] + x[0])/2.0
            xnew = (x[n] + x[0])/2.0
            print "shrinkage"
    x[n] = xnew
    # resort the array
    fx = np.array([f(x_) for x_ in x])
    fxsort = fx.argsort()
    fx = fx[fxsort]
    x = x[fxsort]
    # plot the simplex
    p = plt.Polygon(x, closed=True, fill=False)
    ax = plt.gca()
    ax.add_patch(p)
    fbar = np.sum(fx)/float(n)
    diff_array = np.array([fx[i]-fbar for i in range(0,n+1)])
    standard_error = 1/float(n)*np.sum(diff_array)**2
    if standard_error < epsilon:
        break

bestX = x[0]
bestF = fx[0]

print "number of iterations:"
print count
print "Best variable result:"
print bestX
print "Best function value:"
print bestF

annotation = "f(%3.8f,%3.8f)=%3.8f" % (bestX[0],bestX[1],bestF)
# plot
#plt.axis([-70,150,-150,40])
plt.title("Nelder-Mead Simple 2D Example with f = x^2 + y^2")
plt.xlabel("x")
plt.ylabel("y")
plt.plot(bestX[0],bestX[1],'^')
plt.annotate(annotation, (bestX[0],bestX[1]))
ax2 = plt.gca()
ax2.autoscale_view(True,True,True)
plt.show()
