There are 3 main pieces to our system:
1. The BizHawk emulator. This runs the game that is being played.
2. The Code/ directory. This contains our script as well as the DP1.state file, which is necessary for the system to run.
3. The ROM/ directory. This contains the ROM file, or the game file that will be loaded into BizHawk.

There are some additional and optional pieces:
The Test Files directory contains data from our testing phase, but it is not used to run the system.
The Documentation directory contains our User’s and System manuals.

Startup Procedures:
1. From the BizHawk-1.11.4 directory, open the EmuHawk.exe program.
2. In EmuHawk, select “Open ROM” under the File menu. Then navigate to the “BWDH-FinalProject/ROM” directory and select the “Super Mario World (USA).sfc” file.
3. In EmuHawk, select “Lua Console” under the Tools menu. This will bring up a separate window. From this window select “Open Script…” under the Script menu. Navigate to the “BWDH-FinalProject/Code” directory and select the BWDHmaster.lua file. The project is now running.
