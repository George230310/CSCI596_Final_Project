CC = g++
CFLAGS = -c -I./glm -I./png++ `libpng-config --cflags`
LDFLAGS = -lGL -lGLU -lglut -I./glm -I./png++ `libpng-config --ldflags`

all: raytracer

raytracer: raytracer.o Sphere.o Triangle.o
	$(CC) -o raytracer raytracer.o Sphere.o Triangle.o $(LDFLAGS)

raytracer.o: raytracer.cpp raytracer.h Renderable.h
	$(CC) $(CFLAGS) raytracer.cpp

Sphere.o: Sphere.cpp Renderable.h
	$(CC) $(CFLAGS) Sphere.cpp

Triangle.o: Triangle.cpp Renderable.h
	$(CC) $(CFLAGS) Triangle.cpp

clean:
	rm -f raytracer raytracer.o Sphere.o Triangle.o
