import math
import cmath
import random




def sigmoid(inp):
    # also called logistic curve
    # returns an output between 0 and 1
    # where 0 returns 0.5
    return 1/(1+math.exp(-inp))


def sigmoidBias(inp, bias):
    return sigmoid(inp + bias)


def test():
    assert sigmoid(0.0) == 0.5
    print(sigmoid(1.0))
    print(sigmoid(-1.0))

if __name__=="__main__":
    test()