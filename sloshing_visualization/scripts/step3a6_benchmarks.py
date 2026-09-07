#!/usr/bin/env python3
"""Bounded STEP3A6 stages only; never pilot, full series or CHNS transfer."""
import argparse
from sloshing.multiphase.benchmarks import accuracy_ch

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--stage",required=True,choices=accuracy_ch.STAGES)
    args=parser.parse_args()
    accuracy_ch.STAGES[args.stage]()
