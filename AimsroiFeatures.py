#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 13:44:33 2022

@author: user7t
"""

singularity run -B $HOME/brainvisa-5.0.4:/casa/setup $HOME/brainvisa-5.0.4/brainvisa-5.0.4.sif

export PATH="$HOME/brainvisa-5.0.4/bin:$PATH"

AimsRoiFeatures i-