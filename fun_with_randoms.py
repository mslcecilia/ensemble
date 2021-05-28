# ENGSCI263: Forecasting Geothermal Impacts
# fun_with_randoms.py

# PURPOSE:
# To INTRODUCE you to random numbers.

# PREPARATION:
# TASKS 1 and 2 in main.py

# INSTRUCTIONS:
# turn on/off the code blocks below by toggling True/False 

# import modules and functions
import numpy as np

# Part I: generate and print a random number
    # Run this code several times.
    # i. Does 'a' have the same value each time?
    # ii. What values can 'a' possibly take on?
if True:
    a = np.random.rand()
    print('a =',a)


# Part II: generate repeatable random numbers
    # Run this code several times.
    # i. Does 'a' have the same value each time?
    # ii. Why are 'a' and 'c' the same, but 'b' is different?
if False:
    np.random.seed(1)
    a = np.random.rand()
    print('a =',a)
    b = np.random.rand()
    print('b =',b)
    np.random.seed(1)
    c = np.random.rand()
    print('c =',c)


# Part III: modifying random numbers
    # Run this code several times.
    # i. What values can 'a' possibly take on?
    # ii. What values can 'b' possibly take on?
if False:
    bmin, bmax = [2, 8]
    a = np.random.rand()
    print('a =',a)
    b = a*(bmax-bmin)+bmin
    print('b =',b)


# Part IV: generate normally distributed randoms
    # Run this code several times.
    # i. What values can 'a' possibly take on?
    # ii. What effect does putting an integer inside the randn() function have?
    # iii. What are the statistical properties of randn() numbers?
if False:
    a = np.random.randn()
    print('a =',a)
    a100 = np.random.randn(100)
    print('a100 =',a100)
    print('mean(a100) =', np.mean(a100))
    print('std(a100) =', np.std(a100))

    
# Part V: generate arbitrarily normally distributed randoms
    # Run this code several times.
    # i. What values can 'b100' possibly take on?
    # ii. What are the statistical properties of b100?
if False:
    a100 = np.random.randn(100)
    b100 = a100*5 + 10.
    print('b100 =',b100)
    print('mean(b100) =', np.mean(b100))
    print('std(b100) =', np.std(b100))

