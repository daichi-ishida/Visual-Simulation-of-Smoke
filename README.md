# Visual Simulation of Smoke
This is the implementation of *Visual Simulation of Smoke* by Ronald Fedkiw, Jos Stam, and Henrik Wann Jensen at SIGGRAPH 2001.

## Build
**Windows is not supported.**

Library you need is Eigen 3.3.4

You can build and execute following command.

```shell
$ git clone https://github.com/daichi-ishida/Visual-Simulation-of-Smoke.git
$ make
$ make run
```

When you run this program, simulation result is written to the "output" directory as VTK format.

So to visualize this result, you need **ParaView**.

https://www.paraview.org/

## Screenshot

![Screenshot](output/Screenshot.png)