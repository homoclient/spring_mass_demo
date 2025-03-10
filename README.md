# spring_mass_demo

Needs *Raylib* for visualization, *jo_mepeg.h* for mpeg recording

Everything else in public domain

### 0.c: rope made of 10 particles (1 fixed)

### 1.c, 2.c: solid cube made of particles and springs

press _1_ to rotate, press _2_ to squeeze inward

### 3.c: solid cube falling, colliding with ground

press _space_ to start/pause, press _=_ to start recording, _-_ stop recording. Video is saved as 'rec.mpg', beware the video is vertically flipped

To build on windows:
gcc 3.c raylib.dll -o 3.exe

## Update
### wtf.c: demonstrates a variant of the method which is angular momentum preseving

press _1_ to rotate
