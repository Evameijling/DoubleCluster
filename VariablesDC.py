#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 10:15:16 2021

@author: evagmelichmeijling
"""

DOUBLECLUSTERTXT = open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/AnalysisResults.txt', 'w')
DOUBLECLUSTERTXT.write("VARIABLES:\n")

apertureradius = 10
annulusradius_in = 25
annulusradius_out = 30

DOUBLECLUSTERTXT.write(f"Apertureradius: {apertureradius}\n")
DOUBLECLUSTERTXT.write(f"Annulusradius_in radius: {annulusradius_in}\n")
DOUBLECLUSTERTXT.write(f"Annulusradius_out radius: {annulusradius_out}\n")

DOUBLECLUSTERTXT.write("\n")
DOUBLECLUSTERTXT.close()
