CC = mpicc
CFLAGS = -std=c++11 -Ofast -c -I./libs/glm
LDFLAGS = -lm

all: raytracer

raytracer: raytracer.o Sphere.o Triangle.o
	$(CC) -o raytracer ./libs/lodepng/lodepng.cpp raytracer.o Sphere.o Triangle.o $(LDFLAGS)

raytracer.o: raytracer.cpp raytracer.h Renderable.h
	$(CC) $(CFLAGS) raytracer.cpp

Sphere.o: Sphere.cpp Renderable.h
	$(CC) $(CFLAGS) Sphere.cpp

Triangle.o: Triangle.cpp Renderable.h
	$(CC) $(CFLAGS) Triangle.cpp

clean:
	rm -f raytracer raytracer.o Sphere.o Triangle.o
