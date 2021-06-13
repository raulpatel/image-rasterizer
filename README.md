# image-rasterizer
A software-based computer graphics system that renders imagery via rasterization, including Phong shading, hidden surface remove, and arbitrary camera positions all from the CPU.

# Instructions

1) Install CMake

2) Download VTK and build the files in a directory

3) In CMakeLists.txt, replace the ```<vtk build directory path>``` with the path to your VTK build.

4) From your project directory, enter the command ```cmake .``` to create the Makefile for the project.

5) Type ```make``` to create the executable.


As it is currently set, running the executable will create a directory called ```frames``` with the first frame of the image in rotation. You can change this to whatever number you would like and use a program like ffmpeg to create a movie from the frames.

You can find my example movie here: https://youtu.be/aKMaKMug-ZU
