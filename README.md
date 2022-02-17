### Installation
```
git clone https://github.com/xzh19980906/DMRates.git
cd DMRates
pip install -e .
```

### Description
This is a package to calculate WIMP-xenon scattering event rate, including SI, SD, pion, WIMP-nucleon EFT, WIMP-quark EFT channels.

- CONST.py includes the constants needed in the rate calculation.

- STANDARD.py is the main module to calculate the standard SI or SD event rate. An additional rate for WIMP-pion is included.

- EFT.py is an equivalent copy of MMA WIMP-nucleon EFT package [dmformfactor](https://www.ocf.berkeley.edu/~nanand/software/dmformfactor/). If you have any other better idea to run this package in python, please contact me. See more details in Haxton's paper [arxiv:1308.6288](http://inspirehep.net/record/1251560).

- MATCH.py is a short version of [directdm](https://github.com/DirectDM/directdm-py), which could transfer WIMP-quark EFT operator basis to WIMP-nucleon EFT basis. A next-to-leading order correction to nuclear form factor is added, but still need to be verified. See more details in Zupan+'s paper (can be found [here](https://directdm.github.io/))

- QCD.py gives the parameters needed in MATCH.py.

I offered a few examples for you to get the idea how to use this package.

In the folder Compare_with_EFT, I put the latest comparison between standard spin-dependent theory and effective field theory.

### Attention
Check the parameters in CONST.py before trying to compare the results with other theories, especially halo parameters and dark matter local density etc.

### Contact
If you have any problems, feel free to contact me (xzh19980906@gmail.com).
