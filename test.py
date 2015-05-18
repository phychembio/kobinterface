import math

def stats(x):
  n = 0
  S = 0.0
  m = 0.0
  for x_i in x:
    n = n + 1
    m_prev = m
    m = m + (x_i - m) / n
    S = S + (x_i - m) * (x_i - m_prev)
  return {'mean': m, 'variance': S/n}

def naive_stats(x):
  S1 = sum(x)
  n = len(x)
  S2 = sum([x_i**2 for x_i in x])
  return {'mean': S1/n, 'variance': (S2/n - (S1/n)**2) }

x1 = [1,-1,2,3,0,4.02,5] 
x2 = [x+1e9 for x in x1]

print "naive_stats:"
print naive_stats(x1)
print naive_stats(x2)

print "stats:"
print stats(x1)
print stats(x2)