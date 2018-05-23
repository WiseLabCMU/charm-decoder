# charm-decoder
Charm coherent combining decoder for LoRa packets

# Making sense of existing code

- `finalCountdown10.m`: Use for SF 10
- `finalCountdown.m`: Use for SF 7

Input is the I/Q samples from SX1257
Output is the set of individual SNRs and coherent combined SNR.
Change the file name in line 25 
It is super slow(5 mins) so be patient. 
This also contains lots

- `mycreatechirpspecial.m`: Creates upchirp + downchirp of SF7
- `mycreatechirpspecial1.m`: Creates upchirp + downchirp of SF10/SF12
- `calcSR.m`: Gives signal power