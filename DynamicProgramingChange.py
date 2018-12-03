#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 14:23:48 2018

@author: jennifermichaud
"""

import math
import numpy as np
import copy

def RecursiveChange(money, coins):
    """Using recursion inputs a integer total money amount (money) and
    and a list of coins (coins) and makes change for money using coins using 
    smallest amount of coins possible.
    Outputs the minimum number of coins required. To note this algorithm is
    not efficient.
    """
    if money == 0:
        return 0
    minnumcoins = math.inf
    for i in range(len(coins)):
        if money >= coins[i]:
            numcoins = RecursiveChange(money - coins[i], coins)
            if numcoins + 1 < minnumcoins:
                minnumcoins =  numcoins + 1
    return minnumcoins

                    
def DPChange(money, coins):
    """Using dynamic programing inputs a integer total money amount (money) and
    and a list of coins (coins) and makes change for money using coins using 
    smallest amount of coins possible.
    Outputs the minimum number of coins required and a list of the coins used 
    to reach that total.
    """
    MinNumCoins = np.array([[int(i) for i in range(money+1)], [math.inf for j in range(money+1)]])
    MinNumCoins[1,0] = 0
    coinlistdict = {} # a dictionary to store coins used to make each coin amount
    coinlistdict[0] = []
    for m in range(1,money+1):
        for coin in coins:
            if m >= coin:
                if MinNumCoins[1, m-coin] + 1 < MinNumCoins[1,m]:
                    MinNumCoins[1,m] = MinNumCoins[1,m-coin] +1
                    if coinlistdict[m-coin] != []:
                        coinlistdict[m] = copy.deepcopy(coinlistdict[m-coin])
                        coinlistdict[m].append(coin)   
                    else:
                         coinlistdict[m] = [coin]
    return int(MinNumCoins[1,money]), coinlistdict[money]





##TEST DATA        
#
#money = 11
#coins = [5, 4, 1]                     
#
#
#money = 40
#coins = [50,25,20,10,5,1]    
#
#money = 16509
#coins = [23,20,12,5,3,1] 
#
#money = 21
#coins = [2,3]    
#
#print(RecursiveChange(money, coins))    
#print(DPChange(money, coins))
